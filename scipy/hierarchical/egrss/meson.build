fegrss_source = [
    'src/dtrmm.f90',
    'src/dtrmv.f90',
    'src/dtrsm.f90',
    'src/dtrsv.f90',
]

fegrss_module = custom_target('_fegrss_module',
    input : fegrss_source,
    output : ['_fegrssmodule.c', '_fegrss-f2pywrappers.f'],
    command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', '_fegrss', '--lower', '--build-dir', '@OUTDIR@'],
)

py3.extension_module('_fegrss',
    fegrss_module,
    link_args: version_link_args,
    dependencies: [lapack_dep, blas_dep, fortranobject_dep],
    install: true,
    subdir: 'scipy/hierarchical/egrss'
)

python_sources = [
    '__init__.py',
    '_basic.py',
    'basic.py',
    'egrss.py',
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/hierarchical/egrss',
)

subdir('tests')