egrss_sources = [
  'src/dpotrf.f',
  'src/dtrmm.f',
  'src/dtrmv.f',
  'src/dtrsm.f',
  'src/dtrsv.f',
  'src/dtrtrs.f',
]

# Building ARPACK yields a ton of rank mismatch (scalar and rank-1) warnings
# that cannot be suppressed with a more specific flag.
_suppress_all_warnings = ff.get_supported_arguments('-w')

egrss_lib = static_library('egrss_lib',
  egrss_sources,
  fortran_args: [fortran_ignore_warnings, _suppress_all_warnings],
  include_directories: ['src'],
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

egrss_module = custom_target('egrss_module',
  output: ['_egrssmodule.c', '_egrss-f2pywrappers.f'],
  input: 'egrss.pyf.src',
  command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_egrss = py3.extension_module('_egrss',
  egrss_module,
  link_with: egrss_lib,
  link_args: version_link_args,
  dependencies: [lapack_dep, blas_dep, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/hierarchical/egrss'
)

python_sources = [
    '__init__.py',
    '_basic.py',
    '_decomp_cholesky.py',
    '_egrpack.py',
    'basic.py',
    'egrss.py',
]

py3.install_sources(
    python_sources,
    subdir: 'scipy/hierarchical/egrss',
)

subdir('tests')