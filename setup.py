pkgname="findr"
pkgnamefull="Fast Inference of Networks from Directed Regulations"
version=[0,1,0]
license="AGPL-3"
url="https://github.com/lingfeiwang/findr-python"
urllib="https://github.com/lingfeiwang/findr"
urldoc="https://github.com/lingfeiwang/findr/blob/master/doc.pdf"
author="Lingfei Wang"
author_email="Lingfei.Wang.github@outlook.com"
ftype_np="f"
ftype_c="c_float"
gtype_np="B"
gtype_c="c_ubyte"



def pkg_setup():
	from setuptools import setup
	setup(name=pkgname,
		version='.'.join(map(str,version)),
		author=author,
		author_email=author_email,
		description=pkgnamefull,
		long_description="This is a python interface of "+pkgname+" library. Binary library of "+pkgname+" is required and can be found at "+urllib+". "+pkgname+' ('+pkgnamefull+') is a library to infer gene pairwise regulatory relations and to reconstruct gene regulation networks. A full documentation is available at '+urldoc+'.',
		url=url,
		download_url=url,
		include_package_data=True,
		install_requires=['numpy'],
		classifiers=['Development Status :: 3 - Alpha',
			'Environment :: Console',
			'Intended Audience :: End Users/Desktop',
			'Intended Audience :: Education',
			'Intended Audience :: Science/Research',
			'Operating System :: OS Independent',
			'Programming Language :: Python :: 2.7',
			'Topic :: Scientific/Engineering :: Bio-Informatics',
			'Topic :: Software Development :: User Interfaces'],
		license=license,
		packages=[pkgname],
	)

pkg_setup()
