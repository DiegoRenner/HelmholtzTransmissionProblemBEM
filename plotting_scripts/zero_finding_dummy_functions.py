import numpy as np
import matplotlib.pyplot as plt

# finding zeros of f(x) = min(cos(x),sin(x))
x = np.linspace(-2*np.pi, 2*np.pi,1000)
y = np.min([np.cos(x),np.sin(x)],0)
dy = np.ndarray(x.size)
ddy = np.ndarray(x.size)
flags = np.cos(x) < np.sin(x)
for i,f in enumerate(flags):
    if f:
        dy[i] = -np.sin(x[i])
        ddy[i] = -np.cos(x[i])
    else:
        dy[i] = np.cos(x[i])
        ddy[i] = -np.sin(x[i])

zeros_found = np.array([-3.14166, -1.57074, 3.14161, 4.71242])
values_at_zeros = np.min([np.cos(zeros_found),np.sin(zeros_found)],0)

fig = plt.figure()
plt.plot(x,y)
plt.plot(x,dy)
plt.plot(x,ddy)
plt.scatter(zeros_found, values_at_zeros)
plt.legend(["f(x) = min(cos(x),sin(x))", "f'(x)", "f''(x)", "minima found"], loc="upper right")
#plt.legend(["f(x) = min(cos(x),sin(x))", "cos(x) < sin(x) ? -sin(x) : cos(x)", "cos(x) < sin(x) ? -cos(x) : -sin(x)", "zeros found"], loc="upper right")
plt.savefig("../figures/zeros_dummy_fct1.png")
plt.show()

# finding zeros of f(x) = min(cos(2*x+2),2*sin(x))
x = np.linspace(-2*np.pi, 2*np.pi,1000)
y = np.min([np.cos(2*x+2),2*np.sin(x)],0)
dy = np.ndarray(x.size)
ddy = np.ndarray(x.size)
flags = np.cos(2*x+2) < 2*np.sin(x)
for i,f in enumerate(flags):
    if f:
        dy[i] = -2*np.sin(2*x[i]+2)
        ddy[i] = -4*np.cos(2*x[i]+2)
    else:
        dy[i] = 2*np.cos(x[i])
        ddy[i] = -2*np.sin(x[i])

zeros_found = np.array([-5.71249, -4.14162, -1.5708, 0.570791, 2.14163, 4.71245])
values_at_zeros = np.min([np.cos(2*zeros_found + 2),2*np.sin(zeros_found)],0)

fig = plt.figure()
plt.plot(x,y)
plt.plot(x,dy)
plt.plot(x,ddy)
plt.scatter(zeros_found, values_at_zeros)
plt.legend(["f(x) = min(cos(2*x+2),2*sin(x))", "f'(x)", "f''(x)", "minima found"], loc="upper right")
plt.savefig("../figures/zeros_dummy_fct2.png")
plt.show()



# finding zeros of f(x) = min(cos(2*x+2),2*sin(x))
x = np.linspace(-2*np.pi, 2*np.pi,1000)
y = np.min([np.cos(2*x+2),2*np.sin(x+0.1)],0)
dy = np.ndarray(x.size)
ddy = np.ndarray(x.size)
flags = np.cos(2*x+2) < 2*np.sin(x+0.1)
for i,f in enumerate(flags):
    if f:
        dy[i] = -2*np.sin(2*x[i]+2)
        ddy[i] = -4*np.cos(2*x[i]+2)
    else:
        dy[i] = 2*np.cos(x[i]+0.1)
        ddy[i] = -2*np.sin(x[i]+0.1)

zeros_found = np.array([-5.71249, -4.14164, -1.67082, 0.570823, 2.14157, 4.61249])
values_at_zeros = np.min([np.cos(2*zeros_found+2),2*np.sin(zeros_found+0.1)],0)

fig = plt.figure()
plt.plot(x,y)
plt.plot(x,dy)
plt.plot(x,ddy)
plt.scatter(zeros_found, values_at_zeros)
plt.legend(["f(x) = min(cos(2*x+2),2*sin(x+0.1))", "f'(x)", "f''(x)", "zeros found"], loc="upper right")
plt.savefig("../figures/zeros_dummy_fct3.png")
plt.show()

