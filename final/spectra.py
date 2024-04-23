"""
spectra.py

Plot the energy spectra for case{number}.py

Author: Devon Gulley

Completion date: 2024-04-18

"""
from pyqg import diagnostic_tools as tools
from bt import *

user_choice = input("Enter the case: ")

if user_choice == "1":
    import case1
elif user_choice == "2":
    import case2
elif user_choice == "3":
    import case3
elif user_choice == "4":
    import case4
elif user_choice == "5":
    import case5
else:
    print("No such case")
    quit()


# Following from the example code in the PyQG documentation
m.run()
energy = m.get_diagnostic('KEspec')
enstrophy = m.get_diagnostic('Ensspec')

kr, energy_iso = tools.calc_ispec(m, energy.squeeze())

ks = np.array([3., 80])
es = 5*ks**-4
plt.loglog(kr, energy_iso)
plt.xlabel('Wave Number')
plt.title('Energy Spectrum')
plt.show()
