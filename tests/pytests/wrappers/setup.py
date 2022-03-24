from distutils.core import setup, Extension
import os


source_dir = os.path.abspath('./../../../src')

if __name__ == "__main__":
    setup(
        name="ADDAPyWrappers",
        version="1.0.0",
        description="Python interface for some ADDA functions",
        ext_modules=[Extension(
            "ADDAWrappers", ["wrappers.c", os.path.join(source_dir, "cmplx.c")],
        )]
    )
