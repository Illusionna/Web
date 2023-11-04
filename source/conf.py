project = 'HOME'
copyright = '2023, Illusionna'
author = 'Illusionna'
release = 'v1'


extensions = ['recommonmark', 'sphinx.ext.imgmath']

# extensions = ['recommonmark', 'sphinx.ext.imgmath', 'sphinx.ext.mathjax']



templates_path = ['_templates']
exclude_patterns = []




html_theme = 'alabaster'
html_static_path = ['_static']

html_logo = './logo.png'
html_show_sourcelink = False

imgmath_image_format = 'svg'
imgmath_font_size = 16



# mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML'