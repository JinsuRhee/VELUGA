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

    functions = {key:[] for key in target_functions}
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
            for group, func_names in target_functions.items():
            	if current_function in func_names:
            		
            		functions[group].append((current_function, current_comments))
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
        file.write("IDL Usage\n")
        file.write("=============\n\n")
        file.write("This section explains how to use the VELUGA library with IDL. \n\n")

        for group, funcs in functions.items():
        	file.write(f"{group}\n")
        	file.write("-" * len(group) + "\n\n")

        	for func_name, comments in funcs:
        		file.write(f"{func_name}\n")
        		file.write("-" * len(func_name) + "\n\n")
        		file.write("".join(f"{comment}\n" for comment in comments) + "\n")
        		file.write(".. code-block:: idl\n\n")
        		file.write("    "+f"{func_name}\n\n")
"""
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
"""
# Run

target_functions = {
	"Reading Catalog Functions": ['r_gal'], 
	"Getting Functions": ['g_part'], 
	"Drawing Functions": ['d_2dmap']
}

generate_rst_for_idl('../../veluga__define.pro', 'usage/idl.rst', target_functions)

