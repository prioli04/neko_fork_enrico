import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("coefs.csv")
df.columns = df.columns.str.strip()

print(df[df["t"] > 10]["CD"].mean())

df_downwash = pd.read_csv("output/downwash.csv")
downwash_ellip = - 0.5 / (3.1415 * 10.0)

plt.figure()
plt.xlabel("Time step [-]")
plt.ylabel("CD [-]")
plt.plot(df["CD"], label="Neko ALM (eps=b/8)")
plt.plot([0, len(df["CD"])], [0.00548, 0.00548], "r--", label=r"Matlab ALM (eps=b/8)")
plt.legend()
plt.title("CD convergence")

plt.figure()
plt.xlabel("y [m]")
plt.ylabel(r"$\alpha_i$ [rad]")
plt.title("Downwash")
plt.plot(df_downwash["Points:1"], df_downwash["Velocity:2"], label="Neko ALM result")
plt.plot([-1.0, 1.0], [downwash_ellip, downwash_ellip], "r--", label="Theoretical value for ellip. wing")
plt.grid()
plt.legend()
plt.show()