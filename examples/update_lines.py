from shutil import move, copymode
from os import remove


def update_parameters_in_file(filepath, parameters: dict = None):
    """
        Helper function that updates file parameters to run examples in snakemake
    Args:
        filepath:
        parameters: a dictionary formatted as: {"what to find in line": "what to replace the entire line with"}

    Returns:

    """
    temp = "temp.txt"
    with open(temp, 'w') as new_file:
        with open(filepath) as old_file:
            for line in old_file:
                replaced = False
                for key in list(parameters.keys()):
                    if key in line:
                        new_file.write(parameters[key])
                        replaced = True
                if not replaced:
                    new_file.write(line)
    # Copy the file permissions from the old file to the new file
    copymode(filepath, temp)
    # Remove original file
    remove(filepath)
    # Move new file
    move(temp, filepath)
