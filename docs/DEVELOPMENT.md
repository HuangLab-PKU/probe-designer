# Development Guide

## Contributing

### Commit Message Format

We use [Conventional Commits](https://www.conventionalcommits.org/) for commit messages:

```
<type>[optional scope]: <description>

[optional body]

[optional footer(s)]
```

#### Types:
- `feat`: A new feature
- `fix`: A bug fix
- `docs`: Documentation only changes
- `style`: Changes that do not affect the meaning of the code
- `refactor`: A code change that neither fixes a bug nor adds a feature
- `perf`: A code change that improves performance
- `test`: Adding missing tests or correcting existing tests
- `chore`: Changes to the build process or auxiliary tools

#### Examples:
```bash
git commit -m "feat(probe): add per-gene thermodynamic plots"
git commit -m "fix(mutation): resolve coordinate conversion issues"
git commit -m "docs: update revision records with latest changes"
```

### Branch Naming

- `feature/description`: New features
- `bugfix/description`: Bug fixes
- `hotfix/description`: Critical bug fixes
- `refactor/description`: Code refactoring
- `docs/description`: Documentation updates

### Pull Request Process

1. Create a feature branch from `main`
2. Make your changes with descriptive commits
3. Update documentation if needed
4. Create a pull request with:
   - Clear title and description
   - Reference to any related issues
   - Screenshots for UI changes
   - Test results for functional changes

### Code Review Checklist

- [ ] Code follows project style guidelines
- [ ] Tests pass
- [ ] Documentation is updated
- [ ] No breaking changes (or properly documented)
- [ ] Performance impact considered

## Release Process

### Version Numbering

We follow [Semantic Versioning](https://semver.org/):
- `MAJOR`: Incompatible API changes
- `MINOR`: New functionality in a backwards compatible manner
- `PATCH`: Backwards compatible bug fixes

### Release Steps

1. Update version numbers in relevant files
2. Update CHANGELOG.md with new version
3. Create a git tag: `git tag -a v2.1.0 -m "Release version 2.1.0"`
4. Push tag: `git push origin v2.1.0`
5. Create GitHub release with detailed notes

## Testing

### Running Tests

```bash
# Run all tests
python -m pytest

# Run specific test file
python -m pytest tests/test_mutation_probe.py

# Run with coverage
python -m pytest --cov=src tests/
```

### Test Structure

```
tests/
├── unit/           # Unit tests for individual functions
├── integration/    # Integration tests for workflows
├── fixtures/       # Test data and fixtures
└── conftest.py     # Pytest configuration
```

## Documentation

### Updating Documentation

1. **Code Changes**: Update docstrings and inline comments
2. **API Changes**: Update README.md and API documentation
3. **New Features**: Add to user guide and examples
4. **Bug Fixes**: Update troubleshooting section

### Documentation Structure

- `README.md`: Project overview and quick start
- `CHANGELOG.md`: Detailed change history
- `DEVELOPMENT.md`: This file - development guidelines
- `docs/`: Detailed documentation
- `revision_records.md`: Technical change log

## Environment Setup

### Prerequisites

- Python 3.8+
- Conda or virtual environment
- Git

### Platform Considerations
⚠️ **Important**: RNA structure prediction using ViennaRNA is only available on **Linux** and **macOS** systems. Windows developers should:

1. **Use WSL**: Recommended for Windows development
2. **Use Docker**: For containerized development
3. **Disable RNA structure features**: For testing without RNA structure prediction

### Setup Steps

```bash
# Clone repository
git clone <repository-url>
cd probe_design

# Create conda environment (recommended)
conda env create -f code/environment.yml
conda activate probe_design

# Or create virtual environment with pip
python -m venv probe_design_env
source probe_design_env/bin/activate  # Linux/macOS
# probe_design_env\Scripts\activate   # Windows
pip install -r code/requirements.txt

# For RNA structure prediction (Linux/macOS only)
conda install -c bioconda viennarna

# Install development dependencies
pip install -r requirements-dev.txt
```

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure conda environment is activated
2. **Memory Issues**: Reduce batch size or use smaller test datasets
3. **Performance Issues**: Check thermal filter parameters and sequence lengths

### Getting Help

- Check existing issues on GitHub
- Review documentation
- Create a new issue with detailed description and error logs
