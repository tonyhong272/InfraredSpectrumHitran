from cx_Freeze import setup, Executable
setup( name = "XH",
      version = "1.1",
      description = "plot spectra",
      executables = [Executable("plotspectra.py")]
      )