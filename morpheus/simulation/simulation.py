from morpheus.simulation.cache import SimulationCache
from morpheus.simulation.options import DEFAULT_SIMULATION_OPTIONS, SimulationOptions
from morpheus.simulation.instance import SimulationInstance
from morpheus.interfaces.delta_g import IDeltaG

class Simulation:
    options: SimulationOptions
    cache: SimulationCache

    def __init__(self, options: SimulationOptions=DEFAULT_SIMULATION_OPTIONS):
        self.options = options
        self.cache = SimulationCache()

    def calculate_delta_g(self, obj: IDeltaG, enthalpy = False) -> float:
        instance = SimulationInstance(self.options, self.cache)
        if enthalpy == False:
            instance.result = obj.calculate_delta_g(instance)
            return instance.result
        if enthalpy == True:
            instance.result = obj.calculate_delta_g(instance), obj.extract_delta_h(instance)
            return instance.result
