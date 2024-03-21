import numpy as np
import skfuzzy as fuzz
import matplotlib.pyplot as plt
import pandas as pd

from skfuzzy import control as ctrl
from numpy import genfromtxt

df = pd.read_csv("kidney_diseaseFinalFinal.csv")

bp = ctrl.Antecedent(np.arange(0, 201, 1), 'bp')
bgr = ctrl.Antecedent(np.arange(50, 501, 1), 'bgr')
sc = ctrl.Antecedent(np.arange(0, 31, 1), 'sc')
bu = ctrl.Antecedent(np.arange(0, 501, 1), 'bu')
dm = ctrl.Antecedent(np.arange(0, 1.1, 0.1), 'dm')

ckd = ctrl.Consequent(np.arange(0, 101, 1), 'ckd')



bp['low'] = fuzz.trapmf(bp.universe, [0, 0, 90, 105])
bp['normal'] = fuzz.trimf(bp.universe, [90, 105, 120])
bp['high'] = fuzz.trapmf(bp.universe, [105, 120, 200, 201])
# bp.view()

bgr['safe'] = fuzz.trapmf(bgr.universe, [50, 50, 75, 100])
bgr['borderline'] = fuzz.trapmf(bgr.universe, [75, 100, 125, 150])
bgr['high'] = fuzz.trapmf(bgr.universe, [125, 150, 500, 501])
# bgr.view()

sc['normal'] = fuzz.trapmf(sc.universe, [0, 0, 1.2, 1.8])
sc['high'] = fuzz.trapmf(sc.universe, [1.2, 1.8, 30, 31])
# sc.view()

bu['normal'] = fuzz.trapmf(bu.universe, [0, 0, 45, 55])
bu['high'] = fuzz.trapmf(bu.universe, [45, 55, 500, 501])
# bu.view()

dm['no'] = fuzz.trimf(dm.universe, [0, 0, 1])
dm['yes'] = fuzz.trimf(dm.universe, [0, 1, 1])
# dm.view()

ckd['no'] = fuzz.trapmf(ckd.universe, [0, 0, 45, 55])
ckd['yes'] = fuzz.trapmf(ckd.universe, [45, 55, 100, 101])
# ckd.view()

#2*2*2*3*3 = 72

#RULE non-diabetic patient
rule1 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['no'])
rule2 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['no'])
rule3 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['no'])
rule4 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['no'])

rule5 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['no'])
rule6 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['no'])
rule7 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['high'], ckd['no'])
rule8 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['no'])

rule9 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['high'], ckd['yes'])


#RULE non-diabetic, 1 high
rule10 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['no'])
rule11 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['no'])
rule12 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['no'])
rule13 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['no'])

rule14 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['yes'])
rule15 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['yes'])
rule16 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['low'] & bgr['high'], ckd['yes'])
rule17 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['yes'])

rule18 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['high'] & bgr['high'], ckd['yes'])


rule19 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['low'] & bgr['safe'], ckd['no'])
rule20 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['no'])
rule21 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['no'])
rule22 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['no'])

rule23 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['high'] & bgr['safe'], ckd['yes'])
rule24 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['yes'])
rule25 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['low'] & bgr['high'], ckd['yes'])
rule26 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['high'], ckd['yes'])

rule27 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['high'] & bgr['high'], ckd['yes'])

#RULE non-diabetic, 2 high

rule28 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['low'] & bgr['safe'], ckd['yes'])
rule29 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['yes'])
rule30 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['yes'])
rule31 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['yes'])


#RULE non-diabetic, 3 high

rule32 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['high'] & bgr['safe'], ckd['yes'])
rule33 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['yes'])
rule34 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['low'] & bgr['high'], ckd['yes'])
rule35 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['normal'] & bgr['high'], ckd['yes'])

#RULE non-diabetic, all high

rule36 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['high'] & bgr['high'], ckd['yes'])


#RULE diabetic, low and normal
rule37 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['yes'])
rule38 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['yes'])
rule39 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['yes'])
rule40 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['yes'])

rule41 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['yes'])
rule42 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['yes'])
rule43 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['high'], ckd['yes'])
rule44 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['yes'])

rule45 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['high'], ckd['yes'])


#RULE diabetic, 1 high
rule46 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['yes'])
rule47 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['yes'])
rule48 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['yes'])
rule49 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['yes'])

rule50 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['yes'])
rule51 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['yes'])
rule52 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['low'] & bgr['high'], ckd['yes'])
rule53 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['yes'])

rule54 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['high'] & bgr['high'], ckd['yes'])


rule55 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['low'] & bgr['safe'], ckd['yes'])
rule56 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['yes'])
rule57 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['yes'])
rule58 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['yes'])

rule59 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['high'] & bgr['safe'], ckd['yes'])
rule60 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['yes'])
rule61 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['low'] & bgr['high'], ckd['yes'])
rule62 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['high'], ckd['yes'])

rule63 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['high'] & bgr['high'], ckd['yes'])

#RULE diabetic, 2 high

rule64 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['low'] & bgr['safe'], ckd['yes'])
rule65 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['yes'])
rule66 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['yes'])
rule67 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['yes'])


#RULE diabetic, 3 high

rule68 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['high'] & bgr['safe'], ckd['yes'])
rule69 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['yes'])
rule70 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['low'] & bgr['high'], ckd['yes'])
rule71 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['normal'] & bgr['high'], ckd['yes'])

#RULE diabetic, all high
rule72 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['high'] & bgr['high'], ckd['yes'])



kidney_ctrl = ctrl.ControlSystem([rule1,
                                  rule2, rule3, rule4, rule5, rule6, rule7, rule8, rule9,
                                  rule10, rule11,rule12,
                                  rule13, rule14, rule15, rule16, rule17, rule18, rule19, rule20, rule21, rule22, rule23, rule24,
                                  rule25, rule26, rule27, rule28, rule29, rule30, rule31, rule32, rule33, rule34, rule35, rule36,
                                  rule37, rule38, rule39, rule40, rule41, rule42, rule43, rule44, rule45,
                                rule46, rule47, rule48,
                                  rule49, rule50, rule51, rule52, rule53, rule54, rule55, rule56, rule57, rule58, rule59, rule60,
                                  rule61, rule62, rule63, rule64, rule65, rule66, rule67, rule68, rule69, rule70, rule71,
                                  rule72])

kidney_prediction = ctrl.ControlSystemSimulation(kidney_ctrl)

result_list = []
likelihood_list = []
result_ans = []
predict_diabetes = []

#Iterate each row in the DataFrame

counter = 0
for index, row in df.iterrows():

    kidney_prediction.input['bp'] = row['bp']
    kidney_prediction.input['bgr'] = row['bgr']
    kidney_prediction.input['sc'] = row['sc']
    kidney_prediction.input['bu'] = row['bu']

    if row['dm'] == "yes":
        kidney_prediction.input['dm'] = 1
        predict_diabetes.append(1)
    elif row['dm'] == "no":
        kidney_prediction.input['dm'] = 0
        predict_diabetes.append(0)

    #Evaluate the control system
    kidney_prediction.compute()

    #Obtain output
    likelihood_result = kidney_prediction.output['ckd']
    # print(likelihood_result)

    likelihood_list.append('{:.4f}'.format(likelihood_result))

    if likelihood_result >= 50:
        result_list.append(1)
        result_ans.append("ckd")
    else:
        result_list.append(0)
        result_ans.append("notckd")

# df['likelihood_result'] = result_list
df['result'] = likelihood_list
df['result_classified'] = result_ans
# df['predict_diabetes'] = predict_diabetes

print(df)

# ckd.view(sim=kidney_prediction)

# decoy = input()
comparison = df['classification'] == df['result_classified']

print(comparison.sum())
print(len(df))

similarity_percentage = (comparison.sum() / len(df)) * 100

print(f"Similarity Percentage: {'{:.2f}'.format(similarity_percentage)}%")

# with pd.ExcelWriter('output.xlsx') as excel_writer:
#     df.to_excel(excel_writer, sheet_name='Sheet1', index=False)

# decoy = input()


