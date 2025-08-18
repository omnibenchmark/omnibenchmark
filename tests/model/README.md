# Model Tests

Tests for the new Pydantic-based model that replaces `omni-schema`.

## Structure

- **`test_model.py`**: Core model tests
- **`test_validation_consistency.py`**: Validation behavior tests
- **`test_module.py`**: ModuleMetadata tests (separate concern)
- **`factories.py`**: Factory functions to reduce test boilerplate
- **`conftest.py`**: Common fixtures

## Running

```bash
# All model tests
pytest tests/model -m short

# With coverage
pytest tests/model -m short --cov=omnibenchmark.model
```

## Key Features

- All tests marked with `@pytest.mark.short`
- Factory functions for easy test object creation
- Comprehensive validation coverage
- Legacy API compatibility with deprecation warnings

## Migration Notes

- Replace `omni_schema` imports with `omnibenchmark.model`
- Use direct attribute access instead of getter methods
- Built-in validation via Pydantic