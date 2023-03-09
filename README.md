# Optical Dynamics


A computational physics code for Maxwell-Bloch evolution laser dynamics. 




## Schordinger-Maxwell Equation Dynamics
Inside the optical cavity, the region $0 < z < L-h$ is gain region. We put a passive absorber (with width $h$) and has equi-population-inversion. 

$$\frac{\partial \eta_+(z, t)}{\partial t} = -i{p \over 2\hbar} (n_0E_+ +n_2E_-) -{\eta_+ \over T_2}$$

$$\frac{\partial \eta_-(z, t)}{\partial t} = -i{p \over 2\hbar} (n_0E_- +n_2^{*}E_+) -{\eta_- \over T_2}$$

$$\frac{\partial E_+(z, t)}{\partial t} = {c_0 \over n_r}(iG\eta_+ - lE_+ - {\partial E_+ \over \partial z})$$

$$\frac{\partial E_-(z, t)}{\partial t} = {c_0 \over n_r}(iG\eta_- - lE_- - {\partial E_- \over \partial z})$$

$$\frac{\partial n_0(z, t)}{\partial t} = \frac{n_0 - N_{eq}(z)}{T_1(z)} - \frac{2p}{\hbar} \Im(E_{+}\eta_{+}^{*} +  E_{-} \eta_{-}^{ * } )$$

$$\frac{\partial n_2(z, t)}{\partial t} = -({1 \over T_1} + 4k_w^2 D)n_2 + i {p \over \hbar} (E_+ \eta_{-}^{ * } - E_{-}^{ * } \eta_+)$$




## Partial Differencial Equation Computation Technique

We use reflection correction for boundary condition.

$$E_+(0) = r_l \cdot E_-(0)$$
$$E_-(L) = r_r \cdot E_+(L)$$

and the output:
$$E_{out} = (1 - r_r) \cdot E_+(L)$$  



In the computing part, we mainly use 4th-order Runge-Kutta Algorithm to calculate the dynamic functions for both real and complex part.

$$f( z(t + dt), t + dt) = f(x(t), t) + K dt + o(dt^4)$$

where $K = (k_1 + 2 k_2 + 2  k_3 + k_4)/6$.

$k_1$ is the slope predicted at $t$, i.e. $k_1 = f'(z(t), t)$

$k_2$ is the slope predicted at $t + dt/2$ with $k_1$, i.e. $k_2 = f'(z(t+dt/2), t+dt/2)$

$k_3$ is the slope predicted at $t + dt/2$ with $k_2$, i.e. $k_3 = f'(z(t+dt/2), t+dt/2)$

$k_4$ is the slope predicted at $t + dt$ with $k_3$, i.e. $k_2 = f'(z(t+dt), t+dt)$
