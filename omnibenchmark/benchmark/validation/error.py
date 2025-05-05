class ValidationError(Exception):
    """Exception raised for validation errors."""

    def __init__(self, *errors):
        self.errors = errors[0] if len(errors) == 1 else errors
        super().__init__(*errors)

    def __str__(self):
        if isinstance(self.errors, str):
            return self.errors
        else:
            return "\n".join(map(str, self.errors))
