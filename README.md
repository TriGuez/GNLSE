# Generalised Nonlinear Schrödinger Equation solver for nonlinear optical pulse propagation in optical fibers
This solver is based on the Runge-Kutta 4 in interaction picture algorithm described in [1] for solving the following nonlinear Schrödinger equation [2] : 

$${{\partial A(z,T)} \over {\partial z}} = -{{\alpha(\omega)}\over{2}} A(z,T) +{\sum}_{k \geq 2}{{{i^{k+1}}\over{k!}}\beta_{k}{{\partial^{k}A(z,T)}\over{\partial T^{k}}}} +... $$

$$... i\gamma\Bigg(1+\tau_{shock}{{\partial}\over{\partial T}}\Bigg)\times \Bigg(A(z,T){\int}_{-\infty}^{\infty}{R(T')\vert A(z,T-T')\vert^{2}dT'}\Bigg)\ \ \ \ \ \ \ \ \ \ \textbf{(1)}$$

## Basic usage
Try to run in the Matlab prompt : 
```Matlab
example
```
This script runs the simplest implementation of the solver. Let's study the first section : 
```Matlab
c = 299792458;
Tspan = 4e-12;
lambda_low = 900e-9;
lambda_high = 1300e-9;
l0 = 1040e-9;
f0 = c./l0;
[t, dt, f, df, w, lbd, res] = initGNLSE(Tspan, l0, lambda_low);
wshift = fftshift(w);
tol = 1e-5;
```
This section is creating the usefull axes & values for the simulation, through the `initGNLSE()` function. Also, we need to define some physical values, such as half duration of the time window, the lowest & the highest wavelength of the simulation, and the central frequency of the simulation, which is set here to 1.10 $^{-5}$ .

```Matlab
alpha = 0;
betas = [-3.53571099317077e-26, 3.68095336905838e-41,2.0172694409917e-55...
         ,-1.31835263681886e-69];
gamma = 0.0742400655658342;
L = 1;
fR = 0.18;
```
The next section is defining the optical fibers via 5 parameters : 

* alpha is setting the confinement losses of the fiber, in m $^{-1}$.
* betas is a vector containing the Taylor coefficients of the propagation constant in the fiber, in s $^{n}$.m $^{-1}$ , n $\geq$ 2
* gamma is the nonlinear coefficient of the fiber, in W $^{-1}$.m $^{-1}$. It could also be a vector containing $\gamma(\lambda)$
* L is the fiber length, in m
* fR is the fractionnal Raman response of the fiber. If set to 0, the solver will compute the propagation without considering the Raman response. Typically, $f_R = 0.18$ in fused silica [2].

In the next section, we are defining the input optical pulse : 
```Matlab
N = 2;
tFWHM_p = 100e-15;
t0_p = tFWHM_p/2/sqrt(log(2));
l1 = l0;
f_p = c/l1;
tshift = -3e-12;
P_p = (N^2*abs(betas(1)))./(t0_p^2*gamma);
Epump = sechPulse(P_p,0,t0_p,f_p,tshift,t,f,f0);
```
The pulse is here defined as a 2<sup>nd</sup> order soliton, using the soliton area theorem. $l1$ is defining the central wavelength of the pulse, which is the same here as the central wavelength of the simulation. We also need to compute the corresponding frequency, i.e $f_p$ We can add a temporal shift of the maximum of the pulse. The complex enveloppe of the pulse is then defined using `sechPulse()` as : 

$$ A(z=0,t) = \sqrt{P} \operatorname{sech}\Bigg({{t-t_{shift}}\over{t_0}} \Bigg) \operatorname{exp}\Bigg({i\Big({{C_2}\over{2}}(\omega-\omega_0)^2\Big)-\omega_0 t} \Bigg)\ \ \ \ \ \ \ \ \ \ \textbf{(2)} $$

We can see here that it's possible to add another parameter, called $C_2$ which is corresponding to the second order taylor coefficient of the spectral phase of the pulse. This parameter allow the pulse to be chirped.

Since the physical parameters are defined, we can now focus on the numerical simulation. In the case of a propagation over a small distance, it si possible to create an animation of the propagation. We need first to set some parameters : 

```Matlab
slices = 200;
dL = L/slices;
h = dL/1000;
Esave = zeros(slices+1,length(t));
Eout = Epump;
Esave(1,:) = Epump;
```

The parameters *slices*, *dL* & *h* are used to defined the animation. In this example, we choose to play 200 frames for the propagation over *L*. Then, we need to define the distance to propagate for one frame, which is *dL*. In the end, we need to define the initial stepsize *h* to solve the propagation over *dL*.
Then, we define *Esave* which is a matrix containing the optical field after each propagation over *dL*. The first slice of the matrix is initialized with *Epump*, as well as *Eout* for the solving loop : 
```Matlab
for jk = 2:slices
    Eout = propagationFibre(Eout, dL, h, l0, l0, tol, t, f, lbd,...
             alpha, betas, gamma, fR, 'Propagation example');
    singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'linear')
    title(['L = ' num2str(jk*dL) ' m'])
    Esave(jk,:) = Eout;
    drawnow
end
```
Here, the function `propagationFibre()` is calling the RK4IP solver with the previously defined parameters. It also take as input a string to prompt while computation. The function `singlePlot()` is displaying the resulting optical field in both temporal & spectral domains. At this point, you should get in figure 1 : 

![SSFS](https://user-images.githubusercontent.com/121152666/212467904-056e3cc5-cebf-482a-9ce6-c3bdef9806cc.gif)


In order to understand what happend, we can also use the following line : 
```Matlab
propagationMap(Esave,t,lbd,lambda_low,lambda_high,L);
```
Which will display all the slices in the matrix *Esave* in both domains : 

![image](https://user-images.githubusercontent.com/121152666/217881255-de63f13e-0c2b-43f5-9c73-384f8d56cf86.png)


In the end, we can of course propagates over *L* without animation stuff : 
```Matlab
Eout = propagationFibre(Epump,L,h,l0,l1,tol,t,f,lbd,alpha,betas,gamma,...
    fR, 'Propagation example');
singlePlot(Eout, t, lbd, lambda_low, lambda_high, 'linear')
```

## Solver description
As mentionned before, this solver is based on the Runge-Kutta 4 algorithm for solving ODEs. To use this algortihm on eq. **(1)**, we need to do some maths : 

First, let's rewrite eq. **(1)** as : 

$$ {{\partial A(z,T)} \over {\partial z}} = \Big(\widehat{L} + \widehat{N}\Big) A(z,T)\ \ \ \ \ \ \ \ \ \  \textbf{(3)}$$

Where $\widehat{L}$ stands for the linear part of eq. **(1)**, i.e group velocity dispersion and linear losses operator and $\widehat{N}$ stands for the nonlinear part.

$$\widehat{L} = -{{\alpha(\omega)}\over{2}} +{\sum}_{k \geq 2}{{{i^{k+1}}\over{k!}}\beta_{k}{{\partial^{k}}\over{\partial T^{k}}}}\ \ \ \ \ \ \ \ \ \  \textbf{(4)}$$

$$\widehat{N} = i\gamma\Bigg(1+\tau_{shock}{{\partial}\over{\partial T}}\Bigg)\times\Bigg({\int}_{-\infty}^{\infty}{R(T')\vert A(z,T-T')\vert^{2}dT'}\Bigg)\ \ \ \ \ \ \ \ \ \ \textbf{(5)}$$

Now let's define : 

$$ A_{I}(z,T) = \operatorname{exp}({{(z-z')}\widehat{L}})A(z,T)\ \ \ \ \ \ \ \ \ \ \textbf{(6)}$$

Applied to eq. **(2)** : 

$$ {{\partial A_{I}(z,T)}\over{\partial z}} = \widehat{N_{I}} A_{I}(z,T) \ \ \ \ \ \ \ \ \ \ \textbf{(7)}$$ 

$$ \widehat{N}_{I} = \operatorname{exp}({{-(z-z')}\widehat{L}})\widehat{N} \operatorname{exp}({{(z-z')}\widehat{L}}) \ \ \ \ \ \ \ \ \ \ \textbf{(8)}$$

So, by setting $z' = z+{h\over 2}$, it's easy to integrate eq. **(7)** using standard integration algorithms.

## List of functions

## - adaptiveSolver.m
Usage : 
```Matlab
[Eout, temp_rad, spec_rad] = adaptiveSolver(E, L, h, alpha, betas, gamma, fR,...
                           hR_w, tau_shock, t, lbd, wshift, tol, FiberName)
```
Description :

This function computes the complex enveloppe of an optical pulse after its propagation in a given optical fiber by solving eq. **(1)**. This solver is based on a RK4IP algorithm, and is using an intelligent adaptative stepsize from [3].

Inputs : 

* E : Complex enveloppe of the input optical pulse
* L : Length of the fiber [m]
* h : Initial stepsize [m]
* alpha : Confinement losses of the fiber [m $^{-1}$]
* betas : Taylor coefficients of the propagation constant [s $^{n}$.m $^{-1}$,  n $\geq$ 2]
* gamma : Nonlinear coefficient of the fiber [W $^{-1}$.m $^{-1}$]
* fR : Fractionnal Raman response. Raman scattering is neglected if set to 0
* hR_w : Raman scattering response of the propagation medium
* tau_shock : Shock time for self-steepening [s]
* lbd : Wavelength vector of the simulation [m]
* wshift : Fourier shift of the simulation angular frequency vector [rad.s $^{-1}$]
* tol : Maximum tolerance for the adaptative stepsize
* FiberName : Name of the fiber

Outputs : 

* Eout : Complex enveloppe of the pulse after the fiber
* temp_rad : Duration evolution of the pulse during the propagation [s]
* spec_rad : Spectral width evolution of the pulse during the propagation [m]

## - autocoTrace.m
Usage : 
```Matlab
AC = autocoTrace(E)
```

Description : 

This function computes the autocorrelation function of an optical pulse.

Inputs : 

* E : Complex enveloppe of the input optical pulse OR Intensity enveloppe of the input optical pulse

Outputs : 

* AC : Normalized intensity profile of the autocorrelation trace.

## - compressPulse.m
Usage : 
```Matlab
[Ecomp, Ltot, Dtot] = compressPulse(E, t, f, l0, lc, betas_init, Linit)
```

Description : 

This function calculates the best compressor parameters to minimize the autocorrelation trace width of the input optical field. Initial propagation stepsize & dispersion coefficients must be well guessed in order to reduce the compuation time and to improve the accuracy.

Inputs : 

* E : Complex enveloppe of the optical pulse to compress
* t : Time vector of the simulation [s]
* f : Frequency vector of the simulation [Hz]
* l0 : Central wavelength of the simulation [m]
* lc : Central wavelength of the pulse to compress [m]
* betas_init : Dispersion coefficients [s $^{n}$.m $^{-1}$,  n $\geq$ 2]
* Linit : Initial stepsize for the linear propagation [m]
  
Outputs : 

* Ecomp : Complex enveloppe of the compressed pulse
* Ltot : Total length through the dispersive elements
* Dtot : Total dispersion coefficients for the compressor [s $^{n}$, n $\geq$ 2]

## - gaussianPulse.m
Usage : 
```Matlab
E = gaussianPulse(P,C2,t0,f1,t_shift,t,f, f0)
```

Description :

This function creates the complex enveloppe vector of a gaussian optical pulse. An amplitude noise of one photon per spectral node is added, based on the following : 

$$A_{noise}(t) = \mathcal{F}^{-1}\Bigg[({T_{max}\hbar \omega})^{1/2} \operatorname{exp}(-i\psi(\omega_n))\Bigg] (t) \ \ \ \ \ \ \ \ \ \ \textbf{(9)}$$ 

where $\psi(\omega_n)$ follows a normal distribution.

Inputs : 

* P : Peak power of the Fouier transform limited pulse [W]
* C2 : Second order coefficient of the Taylor developpment of the spectral phase [s $^{2}$]
* t0 : 1/e half pulse duration [s]
* f1 : Optical center frequency [Hz]
* t_shift : Temporal shift of the pulse [s]
* t : Simulation time vector [s]
* f : Simulation frequency vector [Hz]
* f0 : Simulation central frequency [Hz]

Outputs :

* E : Complex enveloppe of the optical pulse

## - gaussRadius.m
Usage : 
```Matlab
radius = gaussRadius(x, y, type)
```

Description : 

This function evaluate the radius of a gaussian function.

Inputs : 

* x : x-axis datas
* y : y-axis datas
* type : string defining the radius (`'1/e'`, `'1/e2'`, `'FWHM'`)

Outpus : 

* radius : Evaluated radius [x-units]

## - initGNLSE.m
Usage : 
```Matlab
[t, dt, f, df, w, lbd, res] = initGNLSE(Tspan, l0, lambda_low)
```

Description : 

Initialize the global values & vectors for the simulation.

Inputs : 

* Tspan : Half width of the temporal window [s]
* l0 : Central wavelength of the simulation [m]
* lambda_low : Lowest wavelength of the simulation [m]

Outputs : 

* t : Created time vector [s]
* dt : Time resolution [s]
* f : Created frequency vector [Hz]
* df : Frequency resolution [Hz]
* w : Created angular frequency vector [rad.s $^{-1}$]
* lbd : Created wavelength vector [m]
* res : Size of the vectors

## - linearProp.m
Usage : 
```Matlab
TE = linearProp(E, f, l0, lc, betas, L)
```

Description : 

Performs the linear propagation of an optical pulse through a dispersive fiber

Inputs : 

* E : Complex enveloppe of the input optical pulse
* f : Frequency vector of the simulation [Hz]
* l0 : Central wavelength of the simulation [m]
* lc : Central wavelength of the optical pulse [m]
* betas : Dispersion coefficients [s $^{n}$.m $^{-1}$,  n $\geq$ 2]
* L : Length of the fiber [m]

Outputs : 

* TE : Complex enveloppe of the output optical pulse

## - nonlinearStep.m
Usage : 
```Matlab
k = nonlinearStep(E, h, gamma, tau_shock, wshift)
```

Description : 

Nonlinear quarter-step for the RK4IP algorithm without Raman scattering (i.e $f_R = 0$).

Inputs : 

* E : Complex enveloppe of the input optical pulse
* h : propagation stepsize [m]
* gamma : Nonlinear coefficient of the fiber [W $^{-1}$.m $^{-1}$]
* tau_shock : Shock time for self-steepening [s]
* wshift : Fourier shift of the simulation angular frequency vector [rad.s $^{-1}$]

Outpus : 

* k : Nonlinear quarter-step operator

## - nonlinearStepFull.m
Usage : 
```Matlab
k = nonlinearStepFull(E, h, fR, hR_w, gamma, tau_shock, wshift)
```

Description : 

Nonlinear quarter-step for RK4IP algorithm including Raman scattering.

Inputs : 

* E : Complex enveloppe of the input optical pulse
* h : propagation stepsize [m]
* fR : Fractionnal Raman response of the propagation medium (typ. 0.18 in fused silica) []
* hR_w : Raman scattering response of the fiber in frequency domain
* gamma : Nonlinear coefficient of the fiber [W $^{-1}$.m $^{-1}$]
* tau_shock : Shock time for self-steepening [s]
* wshift : Fourier shift of the simulation angular frequency vector [rad.s $^{-1}$]

Outpus : 

* k : Nonlinear quarter-step operator

## - propagationFibre.m
Usage : 
```Matlab
[Eout, temp_rad, spec_rad] = propagationFibre(E, L, h, l0, lc, tol ,t, f, lbd, alpha, betas, gamma, fR, FiberName)
```

Description : 


This function computes the propagation over a specific distance of an optical fiber.


INPUTS : 

* E : Complex enveloppe of the input electrical field
* L : Length of the optical fiber [m]
* h : Propagation initial stepsize [m]
* l0 : Central vacuum wavelength of the simulation [m]
* lc : Central vacuum wavelength of the input optical field [m]
* tol : Relative tolerance of the simulation []
* t : Simulation time vector [s]
* f : Simulation frequency vector [Hz]
* lbd : Simulation wavelength vector [m]
* alpha : Confinement losses of the fiber [m $^{-1}$]
* betas : Taylor coefficients of the propagation constant at the given central wavelength [s $^{n}$.m $^{-1}$,  n $\geq$ 2]
* gamma : Nonlinear coefficient of the fiber [W $^{-1}$.m $^{-1}$]
* fR : Fractionnal Raman response of the propagation medium (0.18 in fused silica). Set to zero to exclude Raman scattering []
* FiberName : String to display while computing

OUTPUTS :

* Eout : Complex enveloppe of the output optical field
* temp_rad : pulse duration evolution during the propagation [s]
* spec_rad : Pulse spectral width evolution during the propagation [m]

## - propagationMap.m
Usage : 
```Matlab
propagationMap(E, t, lbd, lambda_low, lambda_high, L)
```

Description : 

Display the slices of a saved propagation matrix in both temporal & spectral domains.

Inputs : 

* E : Matrix of complex enveloppes
* t : time vector of the simulation [s]
* lbd : Wavelength vector of the simulation [m]
* lambda_low : Lowest wavelength to display [m]
* lambda_high : Highest wavelength to display [m]
* L : Propagation length

Outputs : 

## - ramanResponseBW.m

Usage : 
```Matlab
hR_w = ramanResponseBW(fR, wshift)
```

Description : 

This function computes the Raman response of silica acording to [4]

Inputs : 

* fR : Fractionnal Raman response (typ. 0.18 in fused silica)
* wshift : Fourier shift of the simulation angular frequency vector [rad.s $^{-1}$]

Outputs : 

* hR_w : Raman scattering response in the frequency domain

## - rectPulse.m
Usage : 
```Matlab
E = rectPulse(P, tFWHM, f1, t, f, f0)
```

Description : 

This function computes the complex enveloppe of a rectangular optical pulse, with an amplitude noise according to eq. **(9)**.

Inputs : 

* P : Peak power of the pulse [W]
* tFWHM : Full width at half maximum duration of the pulse [s]
* f1 : Central frequency of the pulse [Hz]
* t : Simulation time vector [s]
* f : Simulation frequency vector [Hz]
* f0 : Simulation central frequency [Hz]

Outputs : 

* E : Complex enveloppe of the rectangular pulse

## - RK4IP.m

Usage : 
```Matlab
TE = RK4IP(E, h, alpha, betas, gamma, fR, hR_w, tau_shock, lbd, wshift)
```

Description : 

Runge-Kutta 4 in interaction picture algorithm for the intelligent adaptative stepsize solver for the generalised nonlinear Schrödinger equation **(1)**. This function is the main algorithm of this solver.

Inputs : 

* E : Complex enveloppe of the input electrical field
* h : Propagation initial stepsize [m]
* alpha : Confinement losses of the fiber [m $^{-1}$]
* betas : Taylor coefficients of the propagation constant at the given central wavelength [s $^{n}$.m $^{-1}$,  n $\geq$ 2]
* gamma : Nonlinear coefficient of the fiber [W $^{-1}$.m $^{-1}$]
* fR : Fractionnal Raman response of the propagation medium (0.18 in fused silica). Set to zero to exclude Raman scattering []
* hR_w : Raman scattering response in the frequency domain
* tau_shock : Shock time for self-steepening [s]
* lbd : Simulation wavelength vector [m]
* wshift : Fourier shift of the simulation angular frequency vector [rad.s $^{-1}$]

Outputs : 

* TE : Complex enveloppe of the output pulse

## - sechPulse.m
Usage : 
```Matlab
E = sechPulse(P, C2, t0, f1, t_shift, t, f, f0)
```

Description : 

This function computes the complex enveloppe of an hyperbolic secant optical pulse with amplitude noise according to eq. **(9)**.

Inputs : 

* P : Peak power of the Fouier transform limited pulse [W]
* C2 : Second order coefficient of the Taylor developpment of the spectral phase [s $^{2}$]
* t0 : 1/e half pulse duration [s]
* f1 : Optical center frequency [Hz]
* t_shift : Temporal shift of the pulse [s]
* t : Simulation time vector [s]
* f : Simulation frequency vector [Hz]
* f0 : Simulation central frequency [Hz]

Outputs :

* E : Complex enveloppe of the optical pulse

## - silicaLosses.m
Usage : 
```Matlab
alpha = silicaLosses(lbd)
```

Description : 

This function computes the silica optical losses vs wavelength, according to [5].

Inputs : 

* lbd : Wavelength vector [m]

Outputs : 

* alpha : Silica losses [m $^{-1}$]

## - singlePlot.m
Usage : 
```Matlab
singlePlot(E, t, lbd, lambda_low, lambda_high, spectralScale)
```

Description : 

This function display the input optical field in both time & spectral domain.

Inputs : 

* E : Complex electric field vector
* t : Simulation time vector [s]
* lbd : Simulation wavelength vector [m]
* lambda_low : Lowest wavelength to display [m]
* lambda_high : Highest wavelength to display [m]
* spectralScale : `'linear'` for linear scaling & `'log'` for dB scaling in the spectral domain [string]

Outputs : 



## References
* [1] : Stéphane Balac, Arnaud Fernandez, Fabrice Mahé, Florian Méhats & Rozenn Texier-Picard. *The Interaction Picture method for solving the generalized nonlinear Schrödinger equation in optics*, ESAIM : Mathematical Modelling and Numerical Analysis, 50(4) :945-964, July 2016.
* [2] : Agrawal. *Nonlinear Fiber Optics*, Elsevier, 5<sup>th</sup> edition, 2013.
* [3] : Nguyen, D. T. *Modeling and Design Photonics by Examples using Matlab*, IOP publishing, 2021. doi : 10.1088/978-0-7503-2272-0
* [4] : Blow, K. J., & Wood, D. (1989). *Theoretical description of transient stimulated Raman scattering in optical fibers.* IEEE Journal of Quantum Electronics, 25(12), 2665–2673.
* [5] : Sørensen, S. T. *Deep-blue supercontinuum light sources based on tapered photonic crystal fibers*, PhD thesis, DTU Fotonik ,2013.
