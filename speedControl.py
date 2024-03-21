import numpy as np
import skfuzzy as fuzz
from skfuzzy import control as ctrl

temp = ctrl.Antecedent(np.arange(0, 111, 1), 'temperature')
amb = ctrl.Antecedent(np.arange(0, 101, 1), 'cloud')
speed = ctrl.Consequent(np.arange(0, 101, 1), 'speed')


temp['freezing'] = fuzz.trapmf(temp.universe, [0, 0, 30, 50])
temp['cool'] = fuzz.trimf(temp.universe, [30, 50, 70])
temp['warm'] = fuzz.trimf(temp.universe, [50, 70, 90])
temp['hot'] = fuzz.trapmf(temp.universe, [70, 90, 110, 111])

# temp.view()
# decoy = input()

amb['sunny'] = fuzz.trapmf(amb.universe, [0, 0, 20, 40])
amb['partly_cloudy'] = fuzz.trimf(amb.universe, [20, 50, 80])
amb['overcast'] = fuzz.trapmf(amb.universe, [60, 80, 100, 101])

# amb.view()

speed['slow'] = fuzz.trapmf(speed.universe, [0, 0, 25, 75])
speed['fast'] = fuzz.trapmf(speed.universe, [25, 75, 100, 101])

# speed.view()


rule1 = ctrl.Rule(amb['sunny'] | temp['cool'], speed['fast'])
rule2 = ctrl.Rule(amb['sunny'] | temp['hot'], speed['fast'])
rule3 = ctrl.Rule(amb['partly_cloudy'] | temp['warm'], speed['slow'])
rule4 = ctrl.Rule(amb['partly_cloudy'] | temp['cool'], speed['slow'])

speed_ctrl = ctrl.ControlSystem([rule1, rule2, rule3, rule4])
speeding = ctrl.ControlSystemSimulation(speed_ctrl)

speeding.input['temperature'] = 81.0
speeding.input['cloud'] = 23.0


speeding.compute()

print(speeding.output['speed'])
speed.view(sim=speeding)



decoy = input()
