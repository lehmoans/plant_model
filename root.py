import numpy as np
import matplotlib.pyplot as plt

from plant_model.base import channelBase

class Root(channelBase):
    def __init__(self, R, L, mu, density, theta):
        super().__init__(R, L, mu, density, theta)

    def get_resistance(self):
        # Calculate the resistance of the pipe using the Hagen-Poiseuille equation
        self.resistance = (8 * self.mu * self.L) / (np.pi * self.R**4)
        return self.resistance    
    
    def get_min_pressure_drop(self):
        """
        get minimum pressure difference required to overcome the resistance of the pipe and gravity effects
            min_pressure_diff = resistance* Q {=1} + (density * g * L * sin(theta))
            - assumes positive flow in the direction of the pipe (Pin > Pout)
        
        """
        self.min_pressure_diff = self.get_resistance()+(self.density * 9.81 * self.L * np.sin(np.radians(self.theta)))
        return self.min_pressure_diff

    def get_min_flow_rate(self):
        self.min_flow_rate = (np.pi * self.R**4 * self.get_min_pressure()) / (8 * self.mu * self.L)
        return self.min_flow_rate

    def get_flow_rate(self, P_in=0, P_out=0):
        # Calculate the flow rate using the pressure drop and resistance
        press_drop = P_in - P_out
        if abs(press_drop) < 1e-10:
            print("Forward pressure is negative,flow will be reversed")

        self.flow_rate = -press_drop / self.get_resistance()
        return self.flow_rate
        

    
    

