# Add these settings to enable autodoc for all scripts
autosummary_generate = True  # Generate stub pages for all modules
add_module_names = False     # Remove module names from generated docs
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
    'imported-members': False,
}

# Add the scripts directory to path so Sphinx can find them
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # Root directory
sys.path.insert(0, os.path.abspath('../../scripts'))  # Scripts directory
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
]

# For NumPy/Google style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Theme settings
html_theme = 'sphinx_rtd_theme'