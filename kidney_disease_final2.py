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
bgr['high'] = fuzz.trapmf(bgr.universe, [125, 150, 499, 500])
# bgr.view()

sc['normal'] = fuzz.trapmf(sc.universe, [0, 0, 1.2, 1.8])
sc['high'] = fuzz.trapmf(sc.universe, [1.2, 1.8, 30, 31])
# sc.view()

bu['normal'] = fuzz.trapmf(bu.universe, [0, 0, 45, 55])
bu['high'] = fuzz.trapmf(bu.universe, [45, 55, 500, 501])
# bu.view()

dm['no'] = fuzz.trapmf(dm.universe, [0, 0, 0.3, 0.5])
dm['yes'] = fuzz.trapmf(dm.universe, [0.3, 0.5, 1, 1.1])
# dm.view()

ckd['low'] = fuzz.trapmf(ckd.universe, [0, 0, 25, 40])
ckd['normal'] = fuzz.trimf(ckd.universe, [25, 40,  55])
ckd['high'] = fuzz.trimf(ckd.universe, [40, 55, 70])
ckd['veryhigh'] = fuzz.trapmf(ckd.universe, [55, 70, 100, 101])
# ckd.view()

#RULE non-diabetic, low and normal


# rule0 = ctrl.Rule(dm['yes'], ckd['veryhigh'])

rule1 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['low'])
rule2 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['low'])
rule3 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['low'])
rule4 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['low'])

rule5 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['normal'])
rule6 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['normal'])
rule7 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['high'], ckd['normal'])
rule8 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['normal'])

rule9 = ctrl.Rule(dm['no'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['high'], ckd['high'])


#RULE non-diabetic, 1 high
rule10 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['normal'])
rule11 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['normal'])
rule12 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['normal'])
rule13 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['normal'])

rule14 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['normal'])
rule15 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['normal'])
rule16 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['low'] & bgr['high'], ckd['normal'])
rule17 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['normal'])

rule18 = ctrl.Rule(dm['no'] & sc['high'] & bu['normal'] & bp['high'] & bgr['high'], ckd['high'])


rule19 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['low'] & bgr['safe'], ckd['normal'])
rule20 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['normal'])
rule21 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['normal'])
rule22 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['normal'])

rule23 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['high'] & bgr['safe'], ckd['normal'])
rule24 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['normal'])
rule25 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['low'] & bgr['high'], ckd['normal'])
rule26 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['high'], ckd['normal'])

rule27 = ctrl.Rule(dm['no'] & sc['normal'] & bu['high'] & bp['high'] & bgr['high'], ckd['high'])

#RULE non-diabetic, 2 high

rule28 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['low'] & bgr['safe'], ckd['high'])
rule29 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['high'])
rule30 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['high'])
rule31 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['high'])


#RULE non-diabetic, 3 high

rule32 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['high'] & bgr['safe'], ckd['veryhigh'])
rule33 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['veryhigh'])
rule34 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['low'] & bgr['high'], ckd['veryhigh'])
rule35 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['normal'] & bgr['high'], ckd['veryhigh'])

#RULE non-diabetic, all high

rule36 = ctrl.Rule(dm['no'] & sc['high'] & bu['high'] & bp['high'] & bgr['high'], ckd['veryhigh'])


#RULE diabetic, low and normal 
rule37 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['high'])
rule38 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['high'])
rule39 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['high'])
rule40 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['high'])

rule41 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['high'])
rule42 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['high'])
rule43 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['low'] & bgr['high'], ckd['high'])
rule44 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['high'])

rule45 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['normal'] & bp['high'] & bgr['high'], ckd['veryhigh'])


#RULE diabetic, 1 high
rule46 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['low'] & bgr['safe'], ckd['high'])
rule47 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['safe'], ckd['high'])
rule48 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['low'] & bgr['borderline'], ckd['high'])
rule49 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['borderline'], ckd['high'])

rule50 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['high'] & bgr['safe'], ckd['veryhigh'])
rule51 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['high'] & bgr['borderline'], ckd['veryhigh'])
rule52 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['low'] & bgr['high'], ckd['veryhigh'])
rule53 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['normal'] & bgr['high'], ckd['veryhigh'])

rule54 = ctrl.Rule(dm['yes'] & sc['high'] & bu['normal'] & bp['high'] & bgr['high'], ckd['veryhigh'])


rule55 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['low'] & bgr['safe'], ckd['high'])
rule56 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['high'])
rule57 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['high'])
rule58 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['high'])

rule59 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['high'] & bgr['safe'], ckd['veryhigh'])
rule60 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['veryhigh'])
rule61 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['low'] & bgr['high'], ckd['veryhigh'])
rule62 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['normal'] & bgr['high'], ckd['veryhigh'])

rule63 = ctrl.Rule(dm['yes'] & sc['normal'] & bu['high'] & bp['high'] & bgr['high'], ckd['veryhigh'])

#RULE diabetic, 2 high

rule64 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['low'] & bgr['safe'], ckd['veryhigh'])
rule65 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['normal'] & bgr['safe'], ckd['veryhigh'])
rule66 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['low'] & bgr['borderline'], ckd['veryhigh'])
rule67 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['normal'] & bgr['borderline'], ckd['veryhigh'])


#RULE diabetic, 3 high

rule68 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['high'] & bgr['safe'], ckd['veryhigh'])
rule69 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['high'] & bgr['borderline'], ckd['veryhigh'])
rule70 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['low'] & bgr['high'], ckd['veryhigh'])
rule71 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['normal'] & bgr['high'], ckd['veryhigh'])

#RULE diabetic, all high

rule72 = ctrl.Rule(dm['yes'] & sc['high'] & bu['high'] & bp['high'] & bgr['high'], ckd['veryhigh'])




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
    else:
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

df['likelihood_result'] = result_list
df['result'] = likelihood_list
df['result_classified'] = result_ans
# df['predict_diabetes'] = predict_diabetes

print(df)

# ckd.view(sim=kidney_prediction)

# decoy = input()
comparison = df['classification'] == df['result_classified']

similarity_percentage = (comparison.sum() / len(df)) * 100

print(f"Similarity Percentage: {'{:.2f}'.format(similarity_percentage)}%")

# with pd.ExcelWriter('output_final.xlsx') as excel_writer:
#     df.to_excel(excel_writer, sheet_name='Sheet1', index=False)

# decoy = input()


