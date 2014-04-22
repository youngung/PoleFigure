from distutils.core import setup
setup(name='U-PoleFigure',
      version='0.1',
      description='A simple pole figure plotter',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['upf.analysis','upf.euler','upf.pf','upf.rve','upf.xsym'],
      package_dir={'upf':'src/upf'}
)
