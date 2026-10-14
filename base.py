from abc import ABC, abstractmethod
import numpy as np

class channelBase(ABC):
    def __init__(self, R,L,mu, density, theta):
        self.R = R
        self.L = L
        self.mu = mu
        self.density = density
        self.theta = theta

    @abstractmethod
    def get_resistance(self):
        pass

    @abstractmethod
    def velocity_profile(self):
        pass
