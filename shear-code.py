from symfit import parameters, variables, sin, cos, Fit
import numpy as np
import matplotlib.pyplot as plt

def fourier_series(x, f, n=0):
    """
    Returns a symbolic fourier series of order `n`.

    :param n: Order of the fourier series.
    :param x: Independent variable
    :param f: Frequency of the fourier series
    """
    # Make the parameter objects for all the terms
    a0, *cos_a = parameters(','.join(['a{}'.format(i) for i in range(0, n + 1)]))
    sin_b = parameters(','.join(['b{}'.format(i) for i in range(1, n + 1)]))
    # Construct the series
    series = a0 + sum(ai * cos(i * f * x) + bi * sin(i * f * x)
                     for i, (ai, bi) in enumerate(zip(cos_a, sin_b), start=1))
    return series

x, y = variables('x, y')
w, = parameters('w')
model_dict = {y: fourier_series(x, f=w, n=3)}
print(model_dict)




# Make step function data
eps=0.000000000001
Area=134.65724
xdata = np.array([0.00,0.20,0.40,0.60,0.80,1.00,1.20,1.40,1.60,1.80,2.00]) 
ydata = 160.2176*np.array([0.0,0.22,1.17,2.49,3.87,4.94,5.92,6.78,7.06,7.41,7.60])/(Area*(xdata+eps))

xdata2 = np.array([0.0,0.20,0.40,0.60,0.80,1.00,1.20,1.40,1.60,1.80,2.00]) 
ydata2 = 160.2176*np.array([0.0,0.21,1.05,2.14,3.27,4.27,5.15,5.59,5.99,6.34,6.67])/(Area*(xdata2+eps))


# Define a Fit object for this model and data
fit = Fit(model_dict, x=xdata, y=ydata)
fit_result = fit.execute()
fit2 = Fit(model_dict, x=xdata2, y=ydata2)
fit_result2 = fit2.execute()
print(fit_result)
print(fit_result2)


# Define a Fit object for this model and data
fit = Fit(model_dict, x=xdata, y=ydata)
fit_result = fit.execute()
fit2 = Fit(model_dict, x=xdata2, y=ydata2)
fit_result2 = fit2.execute()
print(fit_result)
print(fit_result2)


dataderivative=((fit.model(x=xdata, **fit_result.params).y)[1:])
data1derivative=((fit.model(x=xdata2, **fit_result2.params).y)[1:])

fitxdata = np.arange(0,2.0,0.001)
fitxdata2 = np.arange(0,2.0,0.001)


# Plot the result
fig, ax = plt.subplots()
ax.set_xlabel('Displacement ($\AA$)', fontsize=15)
ax.set_ylim(0, 7.0)
ax.set_ylabel('Shear Stress (GPa)', fontsize=15)
ax.set_xlim(0, 2.1)
ax.tick_params(axis='both', which='major', labelsize=10, colors='black')
ax.tick_params(axis='both', which='minor', labelsize=10, color = 'black')
plt.plot(fitxdata, fit.model(x=fitxdata, **fit_result.params).y, label='Fitted $<110>$', linewidth=2, color='green')
plt.plot(fitxdata2, fit2.model(x=fitxdata2, **fit_result2.params).y, label='Fitted $<011>$', linewidth=2, color='red')
plt.plot(xdata, ydata, 'go', markersize=10, label='data$<110>$')
plt.plot(xdata2, ydata2, 'r>', markersize=10, label='data$<011>$')
legend = plt.legend(loc='best', shadow=False, fontsize='large')
plt.savefig('Al-alpha-shear.png', format='png', dpi=3000)
plt.show()

print("max of data ",dataderivative.max())

print("max of data 1",data1derivative.max())

