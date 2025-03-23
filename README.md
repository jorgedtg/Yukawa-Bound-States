# Yukawa-Bound-States
Numerical study of bound states and critical behavior of the Yukawa potential, using the Numerov and Runge-Kutta 4 methods. 
The algorithm developed varies the coupling constant λ, the screening parameter α, and the angular momentum quantum number l to calculate the bound states. Critical values of α for each λ were analyzed, revealing their dependence on l.

There are 4 codes in this repository:<br />
  NumerovYukawaL0.m - Numerov's Method to find the number of states and energies for l=0<br />
  NumerovYukawaL1.m - Numerov's Method to find the number of states and energies for l>0<br />
  YukawaACs.m - Algorithm to find the critical values of alpha for lambda = 1-5 and l=0<br />
  YukawaACsLM1.m - Algorithm to find the critical values of alpha for lambda = 1 and l>0<br />

Only the file NumerovYukawaL0.m has detailed comments, but all codes work similarly.
  
I had to do different codes for the l=0 and l>0 cases and to find the critical values of alpha for the following reasons:
  1. In the first two codes, I find the energies' number of states and values. I had to separate them into 2 different codes because I did a shooting method. I found that when the energy I used to integrate the wave function is higher or lower than the real value of the energy, the value of the wave function at r=0 behaves in a particular way, so I started with energy of -30 and added to it until the behavior changed, and after this I lowered the magnitude of the change in energy, for example, if I began adding 1, at the moment it changes behavior, I now add 0.1, and then 0.01, and so on. The behavior for l=0 and l>1 are different, and because the code was already too messy, I decided to separate them into two different codes.
  2. In the two last codes, I look for the critical values of alpha. I separated them into two codes for the reason above. I didn't use the same two codes because I had to search these values in a very big range, and the code I did was too slow for this task. So I decided to try a different method, and instead of obtaining each state for each configuration, I just looked for the ground state and used the bisection method to make it a bit faster. However, I found some other issues that made it slower.
