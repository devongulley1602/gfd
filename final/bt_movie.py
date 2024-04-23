"""
bt_movie.py

Run the time evolution of the BTModel for all case{number}.py

Author: Devon Gulley

Completion date: 2024-04-18

"""
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


# Refresh and run the model
picNum = 0
for snapshot in m.run_with_snapshots(tsnapstart=0, tsnapint=50*m.dt):
    plt.clf()
    plt.tight_layout()
    # the average wind speed in kts
    picNum += 1
    U = (m.u[0] + m.U)*1.94
    V = (m.v[0])*1.94
    magnitude = np.sqrt(U**2 + V**2)
    u_avg = round(np.average(magnitude), 1)
    u_max = round(np.max(magnitude), 1)
    plt.title("PV Anomaly")

    infoString = f'Ld: {"{:.1e}".format(Ld)}m  |  beta/f0: {beta/f0}  |   U0: {U0} kts  |  u_avg: {u_avg} kts  |  u_max: {u_max} kts  |   Days: {np.round(m.t/86400, 2)}  | '

    plt.figtext(0., 0.01, infoString, fontsize=8)

    plt.barbs(X[::spacing, ::spacing], Y[::spacing, ::spacing], U[::spacing, ::spacing], V[::spacing,
              ::spacing], length=4, barbcolor='k', flagcolor='r', linewidth=0.4, sizes=dict(emptybarb=0.00))
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.imshow(m.q.squeeze(), origin='lower', extent=[
               0, m.L, 0, m.W], cmap='RdBu_r')
    plt.pause(1/60)

    plt.draw()
    plt.ioff()
