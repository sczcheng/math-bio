import scipy.integrate
import matplotlib.pyplot as plt

#all in units per day
p2 = 0.2   #probability of vaccination for population #2
ε1 = 0.5    #vaccination site "run-in" rate
ε2 = 0.25
β2 = 8e-11  #contact rate of population #2

θ = 0.00000312   #U.S. birth rate
μ = 0.000024  #U.S. natural death rate 2019
β1 = 8.34e-11  #CA contact rate of population #1 the past 8 weeks
δ = 0.000048 #death due to COVID California the past 8 weeks
γ = 0.00325  #California recovery rate the past 8 weeks
α = 0.012  #immunity loss rate
σ = 0.000003   #relapse rate


def f(t, r):

    S1, I1, R1, S2, I2, R2 = r
    N1 = S1 + R1 + I1
    N2 = S2 + R2 + I2
    p1 = max(p2, (I1+I2)/(S1+R1+I1+S2+R2+I2))

    fS1 = (1-p1)*θ*N1 - β1*S1*I1 - μ*S1 + α*R1 - β2*S1*I2 - p1*ε1*S1
    fI1 = β1*S1*I1 + β2*S1*I2 - μ*I1 - δ*I1 - γ*I1 + σ*R1
    fR1 = p1*ε1*S1 + γ*I1 - σ*R1 - μ*R1 - α*R1 + p1*θ*N1

    fS2 = (1-p2)*θ*N2 - β2*S2*I2 - μ*S2 + α*R2 - β1*S2*I1 - p2*ε2*S2
    fI2 = β2*S2*I2 + β1*S2*I1 - μ*I2 - δ*I2 - γ*I2 + σ*R2
    fR2 = ε2*p2*S2 + γ*I2 - σ*R2 - μ*R2 - α*R2 + p2*θ*N2

    return fS1, fI1, fR1, fS2, fI2, fR2

res = scipy.integrate.solve_ivp(f,(0,100),(0.7,0.2,0.1,0.7,0.2,0.1))

color_map = ["#f0b000", "#CE4867", "#1446C5"]

plt.gcf().set_size_inches(12, 4)

plt.subplot(1, 2, 1)
plt.stackplot(res.t, res.y[:3], labels=['S1', 'I1', 'R1'], colors = color_map)
plt.xlabel('time (day)')
plt.legend()

plt.subplot(1, 2, 2)
plt.stackplot(res.t, res.y[3:], labels=['S2', 'I2', 'R2'], colors = color_map)
plt.legend()
plt.xlabel('time (day)')
plt.tight_layout()

plt.show()
