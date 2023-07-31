**Script for critical current**

to dowload type: git clone https://github.com/EmilRizvanov/current_script.git

Equation are taken from [this article](https://www.researchgate.net/publication/230681050_Temperature_Dependence_of_Pair-breaking_Current_in_Superconductors).

Current_calculations.py the main script which gives csv file with I(T) function.

To use it type: python Current_calculations.py [Number (Amount of Matsubara Frequency)  ] [K Amount of points in the grpah] [Name of csv file with I(T) function]

test.csv was obtained with this command: python3 Current_calculations.py 100 10 test

python plot.py [ name.csv ] make a plot from given csv. (try python3 test.csv)

Red line in the plot in current calculated for formula with works for T close to Tc at T = 0.95Tc.

**!!!**
This solver very sensitive to paramaetrs. And now work only if we put Pi = 3.14, if we put Pi=3.141 it crashes.