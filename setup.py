from distutils.core import setup
setup(name='U-PoleFigure',
      version='0.1',
      description='A simple pole figure plotter',
      author='Youngung Jeong',
      author_email='youngung.jeong@gmail.com',
      packages=['UPF',
                'UPF.conv','UPF.euler','UPF.pf','UPF.pf_analysis','UPF.rve','UPF.xsym'
            ],
      package_dir={'UPF':'src/UPF',
                   'UPF.conv':'src/UPF/conv',
                   'UPF.euler':'src/UPF/euler',
                   'UPF.pf':'src/UPF/pf',
                   'UPF.pf_analysis':'src/UPF/pf_analysis',
                   'UPF.rve':'src/UPF/rve',
                   'UPF.xsym':'src/UPF/xsym'
               },
)
