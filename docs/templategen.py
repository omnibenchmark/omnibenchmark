"""
Generate a yaml template with coments from omnibenchmark models.

It's important to note that the output of this script is intended to be used as
a starting point for creating a benchmark configuration file,
but it's only informative, not normative - the pydantic validations are the
ultimate source of truth.

The generated template includes comments that describe the purpose and expected format of each field.

This is a first implementation, to be improved in the future, and mostly to have a self-documentation feature.
In the future, we should rely on a better template generation mechanism that can
handle more complex scenarios and provide more detailed comments. One
possibility is to extract the information from the model's docstrings, or to
rely more on model's metadata (as in, help_messages or similar)
"""

from omnibenchmark.model.benchmark import Benchmark

try:
    import importlib.metadata

    __version__ = importlib.metadata.version("omnibenchmark")
except ImportError:
    try:
        import pkg_resources

        __version__ = pkg_resources.get_distribution("omnibenchmark").version
    except Exception:
        __version__ = "unknown"


def _get_field_default(model_cls, field_name):
    """Extract the default value for a field from a Pydantic model."""
    try:
        if hasattr(model_cls, "model_fields"):
            field_info = model_cls.model_fields.get(field_name)
            if (
                field_info
                and hasattr(field_info, "default")
                and field_info.default is not ...
            ):
                return field_info.default
    except Exception:
        return None


def generate_yaml_template(model_cls, indent=0, schema_defs=None):
    """Recursively generates a commented YAML template from a Pydantic model."""
    schema = model_cls.model_json_schema()
    if schema_defs is None:
        schema_defs = schema.get("$defs", {})
    yaml_lines = []
    indent_str = "  " * indent

    # Example values for some field names
    example_values = {
        "name": "example benchmark",
        "title": "Example Benchmark Title",
        "description": "A more verbose description, for humans reading this",
        "benchmarker": "Your Name",
        "author": "Your Name",
        "version": "1.0.0",
        "url": "https://example.com",
        "email": "you@example.com",
        # TODO: this one is semantically wrong. In outputs, we actually ask for a "path",
        # but internal validation does not expect a path, neither absolute or
        # relative. It's actually treating it as a templatized filename.
        "path": "filename.txt",
        "command": "example_command",
        "bucket_name": "example-bucket",
        "commit": "c0ffee4",
        "commit_hash": "c0ffee4",
        "benchmark_yaml_spec": "0.0",
        "storage_api": "...",
    }

    for name, prop in schema.get("properties", {}).items():
        # --- Add comments for the field ---
        comments = []
        is_required = name in schema.get("required", [])
        prop_type = prop.get("type", "any")
        if "anyOf" in prop:  # For Optional[T] which becomes anyOf: [T, null]
            non_null_types = []
            for t in prop["anyOf"]:
                if t.get("type") and t.get("type") != "null":
                    non_null_types.append(t.get("type"))
                elif "$ref" in t:
                    ref_name = t["$ref"].split("/")[-1]
                    non_null_types.append(ref_name)
            prop_type = non_null_types[0] if non_null_types else "any"
        elif "$ref" in prop:
            prop_type = prop["$ref"].split("/")[-1]

        # Create main comment with required/optional and type
        main_comment = f"{'required' if is_required else 'optional'} ({prop_type})"
        comments.append(main_comment)

        # Handle defaults - check both direct default and model field defaults
        default_value = prop.get("default")
        if default_value is None:
            # Try to get default from the Pydantic model
            default_value = _get_field_default(model_cls, name)

        if default_value is not None and str(default_value) != "PydanticUndefined":
            comments.append(f"default: {default_value}")

        if prop.get("description"):
            comments.append(prop.get("description"))

        # Add enum options if available
        if "enum" in prop:
            enum_values = prop["enum"]
            comments.append(f"options: {', '.join(map(str, enum_values))}")
        elif "$ref" in prop and prop_type in schema_defs:
            nested_schema = schema_defs[prop_type]
            if "enum" in nested_schema:
                enum_values = nested_schema["enum"]
                comments.append(f"options: {', '.join(map(str, enum_values))}")

        yaml_lines.append(f"{indent_str}# {', '.join(comments)}")

        # --- Add the field entry ---
        # Handle nested models (direct $ref or anyOf with $ref)
        ref_name = None
        if "$ref" in prop:
            ref_name = prop["$ref"].split("/")[-1]
        elif "anyOf" in prop:
            # Handle Optional[Model] which becomes anyOf: [Model, null]
            for option in prop["anyOf"]:
                if "$ref" in option:
                    ref_name = option["$ref"].split("/")[-1]
                    break

        if ref_name:
            yaml_lines.append(f"{indent_str}{name}:")
            if ref_name in schema_defs:
                nested_schema = schema_defs[ref_name]
                if "enum" in nested_schema:
                    # It's an enum, use default if available, otherwise first option
                    if (
                        default_value is not None
                        and str(default_value) != "PydanticUndefined"
                    ):
                        yaml_lines[-1] = f"{indent_str}{name}: {default_value}"
                    else:
                        first_option = nested_schema["enum"][0]
                        yaml_lines[-1] = f"{indent_str}{name}: {first_option}"
                else:
                    nested_lines = _generate_nested_template(
                        nested_schema, indent + 1, schema_defs, example_values
                    )
                    yaml_lines.extend(nested_lines)
            else:
                yaml_lines.append(
                    f"{indent_str}  # ... (nested '{ref_name}' structure) ..."
                )

        # Handle lists of models
        elif prop.get("type") == "array" and "$ref" in prop.get("items", {}):
            ref_name = prop["items"]["$ref"].split("/")[-1]
            yaml_lines.append(f"{indent_str}{name}:")
            yaml_lines.append(f"{indent_str}  - # Example item:")
            if ref_name in schema_defs:
                nested_schema = schema_defs[ref_name]
                nested_lines = _generate_nested_template(
                    nested_schema, indent + 2, schema_defs, example_values
                )
                yaml_lines.extend(nested_lines)
            else:
                yaml_lines.append(f"{indent_str}    # ... ('{ref_name}' object) ...")

        # Handle arrays of simple types
        elif prop.get("type") == "array":
            items_type = prop.get("items", {}).get("type", "any")
            yaml_lines.append(f"{indent_str}{name}:")
            yaml_lines.append(f"{indent_str}  - # Example {items_type}")

        # Handle enums
        elif "enum" in prop:
            enum_values = prop["enum"]
            # Use default if available, otherwise first option
            if default_value is not None and str(default_value) != "PydanticUndefined":
                yaml_lines.append(f"{indent_str}{name}: {default_value}")
            else:
                first_option = enum_values[0] if enum_values else "..."
                yaml_lines.append(f"{indent_str}{name}: {first_option}")

        # Handle simple types
        else:
            # Use default, then example value, then fallback
            if default_value is not None and str(default_value) != "PydanticUndefined":
                example_value = default_value
            elif name.lower() in example_values:
                example_value = example_values[name.lower()]
            else:
                # Try to infer from property type
                if prop_type == "string":
                    example_value = f"example_{name}"
                elif prop_type == "integer":
                    example_value = 0
                elif prop_type == "number":
                    example_value = 0.0
                elif prop_type == "boolean":
                    example_value = False
                else:
                    example_value = "..."

            if example_value is None:
                example_value = ""
            yaml_lines.append(f"{indent_str}{name}: {example_value}")

        # Add line break between fields (except for the last one)
        yaml_lines.append("")

    return "\n".join(yaml_lines)


def _generate_nested_template(schema, indent, schema_defs, example_values=None):
    """Generate template lines for a nested schema."""
    if example_values is None:
        example_values = {}
    yaml_lines = []
    indent_str = "  " * indent

    for name, prop in schema.get("properties", {}).items():
        # Add basic comment
        is_required = name in schema.get("required", [])
        prop_type = prop.get("type", "any")

        if "anyOf" in prop:
            non_null_types = []
            for t in prop["anyOf"]:
                if t.get("type") and t.get("type") != "null":
                    non_null_types.append(t.get("type"))
                elif "$ref" in t:
                    ref_name = t["$ref"].split("/")[-1]
                    non_null_types.append(ref_name)
            prop_type = non_null_types[0] if non_null_types else "any"
        elif "$ref" in prop:
            prop_type = prop["$ref"].split("/")[-1]

        comment = f"{'required' if is_required else 'optional'} ({prop_type})"
        if prop.get("description"):
            comment += f", {prop.get('description')}"

        # Add enum options if available
        if "enum" in prop:
            enum_values = prop["enum"]
            comment += f", options: {', '.join(map(str, enum_values))}"
        elif "$ref" in prop and prop_type in schema_defs:
            nested_schema = schema_defs[prop_type]
            if "enum" in nested_schema:
                enum_values = nested_schema["enum"]
                comment += f", options: {', '.join(map(str, enum_values))}"
        yaml_lines.append(f"{indent_str}# {comment}")

        # Add the field (handle direct $ref or anyOf with $ref)
        ref_name = None
        if "$ref" in prop:
            ref_name = prop["$ref"].split("/")[-1]
        elif "anyOf" in prop:
            # Handle Optional[Model] which becomes anyOf: [Model, null]
            for option in prop["anyOf"]:
                if "$ref" in option:
                    ref_name = option["$ref"].split("/")[-1]
                    break

        if ref_name:
            if ref_name in schema_defs and "enum" in schema_defs[ref_name]:
                # It's an enum, show first option as example
                first_option = schema_defs[ref_name]["enum"][0]
                yaml_lines.append(f"{indent_str}{name}: {first_option}")
            elif ref_name in schema_defs:
                yaml_lines.append(f"{indent_str}{name}:")
                nested_lines = _generate_nested_template(
                    schema_defs[ref_name], indent + 1, schema_defs, example_values
                )
                yaml_lines.extend(nested_lines)
            else:
                yaml_lines.append(f"{indent_str}{name}:")
                yaml_lines.append(
                    f"{indent_str}  # ... (nested '{ref_name}' structure) ..."
                )
        elif prop.get("type") == "array" and "$ref" in prop.get("items", {}):
            ref_name = prop["items"]["$ref"].split("/")[-1]
            yaml_lines.append(f"{indent_str}{name}:")
            yaml_lines.append(f"{indent_str}  - # Example item:")
            if ref_name in schema_defs:
                nested_lines = _generate_nested_template(
                    schema_defs[ref_name], indent + 2, schema_defs, example_values
                )
                yaml_lines.extend(nested_lines)
        elif "enum" in prop:
            enum_values = prop["enum"]
            first_option = enum_values[0] if enum_values else "..."
            yaml_lines.append(f"{indent_str}{name}: {first_option}")
        else:
            # Use default or create a meaningful example
            example_value = prop.get("default")
            if example_value is None:
                # Use example values mapping first, then fallback
                if name.lower() in example_values:
                    example_value = example_values[name.lower()]
                elif name.lower() in ["name", "title"]:
                    example_value = f"example_{name}"
                elif name.lower() in ["description"]:
                    example_value = f"Example {name} description"
                elif prop_type == "string":
                    example_value = f"example_{name}"
                elif prop_type == "integer":
                    example_value = 1
                elif prop_type == "number":
                    example_value = 1.0
                elif prop_type == "boolean":
                    example_value = True
                else:
                    example_value = "..."

            if example_value is None:
                example_value = ""
            yaml_lines.append(f"{indent_str}{name}: {example_value}")

        # Add line break between fields (except for the last one)
        yaml_lines.append("")

    return yaml_lines


if __name__ == "__main__":
    import os

    template = generate_yaml_template(Benchmark)

    # Save it inside the templates directory
    output_path = "templates/benchmark_template.yaml"

    # Ensure the templates directory exists
    os.makedirs("templates", exist_ok=True)

    # Write the template to file
    with open(output_path, "w") as f:
        f.write("# Omnibenchmark YAML Template.\n")
        f.write(
            f"# Generated from Pydantic Benchmark model (omnibenchmark version: {__version__})\n"
        )
        f.write("# DO NOT EDIT THIS FILE\n")
        f.write("")
        f.write(template)

    print(f"✅ Generated template at {output_path}")
