from pathlib import Path

import dendropy
import typer


def get_schema(path: Path):
    if path.suffix in [".nxs", ".nexus", ".treefile"]:
        return "nexus"
    if path.suffix in [".newick", ".nwk"]:
        return "newick"
    raise Exception(f"Cannot get schema for file {path}")


def tree_converter(
    input: Path = typer.Argument(..., help="The path to the tree file in newick format."),
    output: Path = typer.Argument(..., help="The path to the tree file in newick format."),
    input_schema: str = typer.Option(
        "", help="The input file schema. If empty then it tries to infer the schema from the file extension."
    ),
    output_schema: str = typer.Option(
        "", help="The output file schema. If empty then it tries to infer the schema from the file extension."
    ),
    node_label: str = typer.Option("", help="Label the nodes from an annotation if present."),
):
    if not input_schema:
        input_schema = get_schema(input)
    if not output_schema:
        output_schema = get_schema(output)

    input_str = input.read_text()

    # Replace quote char because dendropy nexus tokenizer only uses single quotes by default
    input_str = input_str.replace('"', "'")

    tree = dendropy.Tree.get(data=input_str, schema=input_schema)

    if node_label:
        for node in tree:
            try:
                node.label = "%.2g" % float(node.annotations.get_value(node_label))
            except:
                node.label = node.annotations.get_value(node_label)

    tree.write(path=output, schema=output_schema, suppress_rooting=True)


if __name__ == "__main__":
    typer.run(tree_converter)