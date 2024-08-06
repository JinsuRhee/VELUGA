import os

def parse_idl_functions(file_path, target_functions):
    """
    Parse IDL file to extract function or procedure definitions.
    
    Parameters
    ----------
    file_path : str
        Path to the IDL file to parse.
    
    Returns
    -------
    dict
        Dictionary with keys 'r_*' and 'g_*' mapping to lists of function definitions.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    functions = {'r_*': [], 'g_*': [], 'd_*': []}
    current_function = None
    current_comments = []

    for line in lines:
        line = line.strip()
        if line.startswith(';'):
            # Accumulate comment lines
            current_comments.append(line.lstrip(';').strip())
        elif line.startswith('PRO ') or line.startswith('FUNCTION '):
            # Function or procedure definition
            if not '::' in line: continue
            current_function = line.split()[1].split('::')[1]  # Handle functions with arguments
            current_function = current_function[:-1]
            print(current_function)
            if current_function in target_functions:
            	if current_function.startswith('r_'):
            		functions['r_*'].append((current_function, current_comments))
            	elif current_function.startswith('g_'):
                	functions['g_*'].append((current_function, current_comments))
            current_comments = []  # Reset comments for next function

    return functions

def generate_rst_for_idl(file_path, output_rst, target_functions):
    """
    Generate RST documentation for IDL functions.
    
    Parameters
    ----------
    file_path : str
        Path to the IDL file.
    output_rst : str
        Path where the generated RST file will be saved.
    """
    functions = parse_idl_functions(file_path, target_functions)

    with open(output_rst, 'w') as file:
        file.write("IDL Functions\n")
        file.write("=============\n\n")

        file.write("Reading Functions\n")
        file.write("-----------------\n\n")
        for func_name, comments in functions['r_*']:
            file.write(f"{func_name}\n")
            file.write("-" * len(func_name) + "\n\n")
            file.write("".join(f"{comment}\n" for comment in comments) + "\n")
            file.write(".. code-block:: idl\n\n")
            file.write(f"    {func_name}\n\n")

        file.write("Getting Functions\n")
        file.write("-----------------\n\n")
        for func_name, comments in functions['g_*']:
            file.write(f"{func_name}\n")
            file.write("-" * len(func_name) + "\n\n")
            file.write("".join(f"{comment}\n" for comment in comments) + "\n")
            file.write(".. code-block:: idl\n\n")
            file.write(f"    {func_name}\n\n")

        file.write("Drawing Functions\n")
        file.write("-----------------\n\n")
        for func_name, comments in functions['d_*']:
            file.write(f"{func_name}\n")
            file.write("-" * len(func_name) + "\n\n")
            file.write("".join(f"{comment}\n" for comment in comments) + "\n")
            file.write(".. code-block:: idl\n\n")
            file.write(f"    {func_name}\n\n")
# Run

target_functions = [
	'r_gal', 
	'g_part', 
	'd_2dmap'
]

generate_rst_for_idl('../../veluga__define.pro', 'usage/idl.rst', target_functions)

