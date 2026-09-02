"""Gene-aware BLAST submission, parsing, and filter.

Used by mutation (and TCR, next round) pipelines. The filter rule:

    Reject a probe only if it has a 100%-identity hit to a DIFFERENT gene,
    OR three or more >=95% hits to different genes.
    Same-gene hits across isoforms / species are accepted (they reflect
    paralog/ortholog identity, not off-target binding in a human sample).

This is distinct from ``probe_designer.filtering.SequenceFilter`` (which
applies a *species-aware* filter for the mRNA pipeline). The two coexist;
they answer different specificity questions.
"""
from __future__ import annotations

import concurrent.futures
import logging
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

logger = logging.getLogger(__name__)


# Known full-name gene synonyms for the same-gene filter. The list is small;
# extend in-place when new false rejects are observed in lab BLAST output.
GENE_SYNONYMS: Dict[str, List[str]] = {
    "GSN": ["gelsolin"],
    "ACTB": ["actin, beta", "beta-actin", "beta actin", "actin beta"],
}


# ----------------------------------------------------------------------
# Pure helpers
# ----------------------------------------------------------------------


def is_same_gene_hit(subject_def: str, gene_symbol: str) -> bool:
    """Return ``True`` when ``subject_def`` (BLAST hit description) names the
    same gene as ``gene_symbol`` (any isoform, any species).

    Matches via:
      - parenthesised symbol ``(SYM)``,
      - whole-word symbol case-insensitive,
      - registered full-name synonym in :data:`GENE_SYNONYMS`.
    """
    if not subject_def or not gene_symbol:
        return False
    sym = gene_symbol.upper()
    subj_low = subject_def.lower()
    # Case-insensitive: NCBI names orthologs / read-through loci with mixed case
    # (e.g. rodent "Tp53", the human "C12orf57" locus); an exact-case match would
    # miss the same-gene hit and falsely reject the probe as off-target.
    patterns = [rf"\({sym}\)", rf"\b{sym}\b"]
    for pat in patterns:
        if re.search(pat, subject_def, re.IGNORECASE):
            return True
    for syn in GENE_SYNONYMS.get(sym, []):
        if syn in subj_low:
            return True
    return False


def extract_binding_query(arm5: str, arm3: str, chemistry: str) -> str:
    """Build the BLAST query string for one probe.

    For non-iLock chemistries the binding region is the literal arm3+arm5
    concatenation. For iLock the SNP nt is shared between ``arm3[-1]`` and
    ``arm5[0]`` (invader-overlap); concatenating verbatim would feed BLAST a
    chimeric sequence the transcriptome never contains, so the duplicate is
    dropped: query = ``arm3[:-1] + arm5``.
    """
    if not arm3 or not arm5:
        return ""
    if chemistry == "iLock":
        return (arm3[:-1] + arm5).upper()
    return (arm3 + arm5).upper()


def binding_query_junction_index(arm3: str) -> int:
    """0-based index in the binding query where the arm5 side begins.

    The ligation junction sits between ``query[idx-1]`` and ``query[idx]``, so
    ``idx`` is simply the number of bases on the arm3 side. That is
    ``len(arm3)`` for **both** chemistries: for iLock,
    ``extract_binding_query`` drops the duplicated overlap base, and since
    ``arm5[0] == arm3[-1]`` the string ``arm3[:-1] + arm5`` still carries a full
    ``arm3`` in its first ``len(arm3)`` positions.
    """
    return len(arm3)


def hsp_spans_junction(
    query_from: int, query_to: int, junction_index: int, clamp_nt: int,
) -> bool:
    """Whether an HSP covers the ligation junction with the ligase's clamp.

    ``query_from`` / ``query_to`` are BLAST's 1-based inclusive query
    coordinates. The hit must cover ``clamp_nt`` bases on each side of the
    junction — anything less cannot be closed by a ligase, so it is an
    off-target the probe *binds*, not one it can circularise on.

    This distinction is the whole point: a probe whose arm happens to match 15/15
    somewhere is not a false-signal risk, while a 39/40 hit straddling the
    junction very much is. Judging on identity alone conflates them.
    """
    if query_from <= 0 or query_to <= 0:
        return False
    start0 = query_from - 1                 # 0-based inclusive
    end0 = query_to                         # 0-based exclusive
    return start0 <= junction_index - clamp_nt and end0 >= junction_index + clamp_nt


# ----------------------------------------------------------------------
# Batched BLAST submission + caching
# ----------------------------------------------------------------------


def run_batch_blast_search(
    seqs: List[Dict],
    blast_dir: Path,
    label: str,
    *,
    batch_size: int = 100,
    poll_interval_s: int = 360,
    max_polls: int = 10,
    max_workers: int = 5,
) -> bool:
    """Submit ``seqs`` to NCBI BLAST in batches; cache one XML per batch.

    Each ``seqs`` element must carry ``probe_id`` and ``binding_sequence``.
    Resumable: any batch whose XML already exists in ``blast_dir`` is skipped.

    ``label`` is prepended to log messages so callers running multiple
    chemistries / panels in parallel can tell their batches apart.

    Returns ``True`` if at least one batch XML was produced (cached or new).
    """
    blast_dir = Path(blast_dir)
    blast_dir.mkdir(parents=True, exist_ok=True)
    batches = [seqs[i:i + batch_size] for i in range(0, len(seqs), batch_size)]
    logger.info("[%s] %d sequences -> %d batches", label, len(seqs), len(batches))

    def process(batch: List[Dict], idx: int) -> Optional[Path]:
        xml_path = blast_dir / f"batch_{idx + 1}.xml"
        if xml_path.exists():
            logger.info("[%s] Batch %d cached", label, idx + 1)
            return xml_path
        combined = "\n".join(f">{s['probe_id']}\n{s['binding_sequence']}" for s in batch)
        params = {
            "CMD": "Put", "PROGRAM": "blastn", "DATABASE": "refseq_rna",
            "QUERY": combined, "HITLIST_SIZE": "10", "EXPECT": "1e-5",
            "WORD_SIZE": "7", "GAPCOSTS": "5 2", "PENALTY": "-3",
            "REWARD": "2", "FORMAT_TYPE": "XML",
        }
        logger.info("[%s] Submitting batch %d/%d (%d seqs)...",
                    label, idx + 1, len(batches), len(batch))
        r = requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=params)
        r.raise_for_status()
        rid = None
        for line in r.text.split("\n"):
            if "RID =" in line:
                rid = line.split("RID =")[1].strip()
                break
        if not rid:
            logger.error("[%s] No RID for batch %d", label, idx + 1)
            return None
        logger.info("[%s] Batch %d RID: %s", label, idx + 1, rid)

        for attempt in range(max_polls):
            time.sleep(poll_interval_s)
            s = requests.get(
                "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"},
            )
            s.raise_for_status()
            if "<?xml version=" in s.text and "<BlastOutput>" in s.text:
                xml_path.write_text(s.text, encoding="utf-8")
                logger.info("[%s] Batch %d complete", label, idx + 1)
                return xml_path
            logger.info("[%s] Batch %d waiting... (%d/%d)",
                        label, idx + 1, attempt + 1, max_polls)
        logger.error("[%s] Batch %d timed out", label, idx + 1)
        return None

    cached = {i for i in range(len(batches))
              if (blast_dir / f"batch_{i + 1}.xml").exists()}
    missing = [i for i in range(len(batches)) if i not in cached]
    if cached:
        logger.info("[%s] %d batches already cached", label, len(cached))
    if not missing:
        return True

    workers = min(max_workers, len(missing))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(process, batches[i], i): i for i in missing}
        for fut in concurrent.futures.as_completed(futures):
            idx = futures[fut]
            try:
                fut.result()
            except Exception as e:
                logger.error("[%s] Batch %d exception: %s", label, idx + 1, e)

    done = sum(1 for i in range(len(batches))
               if (blast_dir / f"batch_{i + 1}.xml").exists())
    logger.info("[%s] BLAST completed: %d/%d batches", label, done, len(batches))
    return done > 0


# ----------------------------------------------------------------------
# Result parsing
# ----------------------------------------------------------------------


def parse_blast_results(
    seqs: List[Dict],
    blast_dir: Path,
    *,
    batch_size: int = 100,
) -> Dict[str, List[Dict]]:
    """Parse the per-batch XML files written by :func:`run_batch_blast_search`.

    Returns a dict ``{probe_id: [hit_dict, ...]}`` where each hit dict has
    ``subject_id``, ``align_length``, ``identities``, ``match_percentage``,
    ``evalue``.
    """
    blast_dir = Path(blast_dir)
    results: Dict[str, List[Dict]] = {}
    n_batches = (len(seqs) + batch_size - 1) // batch_size
    for bi in range(1, n_batches + 1):
        f = blast_dir / f"batch_{bi}.xml"
        if not f.exists():
            logger.warning("Missing batch %d XML at %s", bi, f)
            continue
        start = (bi - 1) * batch_size
        batch_seqs = seqs[start:start + batch_size]
        try:
            tree = ET.parse(f)
            iters = tree.getroot().findall(".//Iteration")
            for k, it in enumerate(iters):
                if k >= len(batch_seqs):
                    break
                pid = batch_seqs[k]["probe_id"]
                hits: List[Dict] = []
                hit_list_node = it.find("Iteration_hits")
                if hit_list_node is not None:
                    for hit in hit_list_node.findall("Hit"):
                        hsps = hit.find("Hit_hsps")
                        if hsps is None:
                            continue
                        for hsp in hsps.findall("Hsp"):
                            ids_el = hsp.find("Hsp_identity")
                            al_el = hsp.find("Hsp_align-len")
                            if ids_el is None or al_el is None:
                                continue
                            ids = int(ids_el.text)
                            al = int(al_el.text)
                            pct = (ids / al) * 100 if al > 0 else 0.0
                            ev = hsp.find("Hsp_evalue")
                            qf = hsp.find("Hsp_query-from")
                            qt = hsp.find("Hsp_query-to")
                            hits.append({
                                "subject_id": (hit.find("Hit_def").text or ""),
                                "align_length": al,
                                "identities": ids,
                                "match_percentage": pct,
                                "evalue": float(ev.text) if ev is not None else 1.0,
                                # 1-based inclusive, BLAST's own convention. Needed
                                # to ask whether the hit covers the ligation
                                # junction — an off-target the probe merely binds
                                # is a different problem from one it can close on.
                                "query_from": int(qf.text) if qf is not None else 0,
                                "query_to": int(qt.text) if qt is not None else 0,
                            })
                results[pid] = hits
        except Exception as e:
            logger.error("Parse error batch %d: %s", bi, e)
    return results


# ----------------------------------------------------------------------
# Filter + dedup
# ----------------------------------------------------------------------


#: Identity a hit needs before it counts against a probe at all.
DEFAULT_IDENTITY_THRESHOLD: float = 95.0
#: Aligned length a hit needs before it counts. Below roughly this the duplex is
#: not stable at the hybridisation temperature, so a short perfect match inside
#: one arm is noise. The old rule had no length term and rejected a good probe on
#: a 15/15 arm-internal hit.
DEFAULT_MIN_ALIGN_LENGTH: int = 20
#: Contiguous bases the ligase needs on each side of the nick — see
#: ``qc.cross_ligation.DEFAULT_VICINITY_N`` and [[reference-splintr-fidelity]].
DEFAULT_CLAMP_NT: int = 3
#: How many non-junction off-target binders are tolerated before the probe is
#: judged promiscuous. These do not miscall, they sequester probe.
DEFAULT_MAX_OFF_TARGET_BINDING: int = 3


def apply_gene_aware_filter(
    seqs: List[Dict],
    blast: Dict[str, List[Dict]],
    *,
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
    min_align_length: int = DEFAULT_MIN_ALIGN_LENGTH,
    clamp_nt: int = DEFAULT_CLAMP_NT,
    max_off_target_binding: int = DEFAULT_MAX_OFF_TARGET_BINDING,
) -> Tuple[List[Dict], List[Dict]]:
    """Split ``seqs`` into (passed, rejected), judging hits by ligation competence.

    Each element must carry ``probe_id``, ``gene``, and either ``junction_index``
    or ``arm3`` (so the junction can be located in the query).

    **What changed, and why.** The rule used to be "reject on any 100%-identity
    off-target hit, or on three >=95% hits", with no reference to alignment
    length or to *where* on the query the hit fell. That is wrong in both
    directions. A 15/15 perfect hit sitting inside one arm rejected a perfectly
    good probe, because identity was computed over the HSP rather than the
    footprint. Meanwhile a 39/40 hit straddling the ligation junction — an
    off-target the probe will circularise on and report as real signal — scored
    97.5%, not 100%, so it needed two more like it before anything happened.
    That second failure mode is the one eLife 107070 documented in 10x's own
    Xenium panels.

    So a hit is now judged on whether a **ligase could close on it**:

    * ``off_target_ligation_competent`` — spans the junction with ``clamp_nt``
      bases each side, at or above ``identity_threshold``, at least
      ``min_align_length`` long. One is enough to reject: it produces signal
      indistinguishable from the real target.
    * ``off_target_binding`` — same identity and length, but not spanning the
      junction. It cannot miscall, only sequester probe, so it takes
      ``max_off_target_binding`` of them to reject.

    Passed records gain ``blast_hits``, ``best_match_percentage``,
    ``same_gene_exact``, and ``off_target_binding_hits``. Rejected records carry
    ``reason`` and ``off_target_subjects``.

    Raises ``ValueError`` if a record cannot supply its junction index — the
    junction rule silently degrading to the old identity-only behaviour is
    exactly how this bug survived.
    """
    passed: List[Dict] = []
    rejected: List[Dict] = []
    for s in seqs:
        pid = s["probe_id"]
        gene = s["gene"]
        if "junction_index" in s:
            junction = int(s["junction_index"])
        elif s.get("arm3"):
            junction = binding_query_junction_index(s["arm3"])
        else:
            raise ValueError(
                f"probe {pid}: BLAST record carries neither 'junction_index' nor "
                f"'arm3', so the ligation junction cannot be located in the query "
                f"— refusing to fall back to an identity-only verdict"
            )

        alignments = blast.get(pid, [])
        exact = [a for a in alignments if a.get("match_percentage", 0) >= 100.0]
        credible = [
            a for a in alignments
            if a.get("match_percentage", 0) >= identity_threshold
            and a.get("align_length", 0) >= min_align_length
        ]
        off_target = [
            a for a in credible
            if not is_same_gene_hit(a.get("subject_id", ""), gene)
        ]
        ligatable = [
            a for a in off_target
            if hsp_spans_junction(
                a.get("query_from", 0), a.get("query_to", 0), junction, clamp_nt,
            )
        ]
        binding_only = [a for a in off_target if a not in ligatable]

        if ligatable:
            rejected.append({
                **s,
                "reason": "off_target_ligation_competent",
                "off_target_subjects": [a["subject_id"][:60] for a in ligatable[:3]],
            })
        elif len(binding_only) >= max_off_target_binding:
            rejected.append({
                **s,
                "reason": "multiple_off_target_high_similarity",
                "off_target_subjects": [a["subject_id"][:60]
                                        for a in binding_only[:3]],
            })
        else:
            best = max((a.get("match_percentage", 0) for a in alignments), default=0)
            same_gene_exact = sum(
                1 for a in exact if is_same_gene_hit(a.get("subject_id", ""), gene)
            )
            passed.append({
                **s,
                "blast_hits": len(alignments),
                "best_match_percentage": best,
                "same_gene_exact": same_gene_exact,
                "off_target_binding_hits": len(binding_only),
            })
    return passed, rejected


def find_best_per_mutation_record(
    passed: List[Dict],
    *,
    key_fields: Tuple[str, ...] = (
        "Gene.refGene.x", "Start", "End", "Ref", "Alt",
    ),
) -> List[Dict]:
    """Group passed BLAST records by the mutation key and return the lowest-score
    candidate per mutation. Each input record must carry ``original_item``
    (the mutation record), ``original_probe`` (one of its ``plp_mut``), and
    must inherit ``score`` from the probe.

    Returns one dict per mutation: ``{seq_info, item, probe, score}``.
    """
    groups: Dict[str, List[Dict]] = {}
    for s in passed:
        item = s["original_item"]
        probe = s["original_probe"]
        key = "_".join(str(item.get(field, "")) for field in key_fields)
        groups.setdefault(key, []).append({
            "seq_info": s,
            "item": item,
            "probe": probe,
            "score": probe.get("score", float("inf")),
        })
    return [min(g, key=lambda x: x["score"]) for g in groups.values()]
