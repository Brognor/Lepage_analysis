import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['text.usetex'] = True 
energy_levels=[1,2,3,4,5,6,7,8,9,10]
data = np.array([-1.568949087013607, -0.19602247419243213, -0.07505287685489748, -0.039067748657544143, -0.023876894192653708, -0.016084664603113197, -0.011565622189664282, -0.008714093928574584, -0.006800400515203364, -0.0054541801364393905])

coulomb=[-1/(2*n**2) for n in range(1,11)]

delta_1st_order_fit =np.array([-1.5603378734445712, -0.2575422340295953, -0.09482732869675746, -0.04781777925915638, -0.028482702939593457, -0.018797860490510812, -0.013295445704167864, -0.009883472365557791, -0.007627349576642181, -0.0060603000875485875])
delta_1st_order_low =np.array([-0.9541801065643083, -0.18177251316956242, -0.0723770410345255, -0.038346564151652274, -0.023633440804551355, -0.015991574532731815, -0.011528221894021324, -0.008699570477119778, -0.0067958574410039605, -0.005454142320668325])

c_delta_fit_a1 =np.array([-1.5689291210037482, -0.1977465403797396, -0.07344981031565112, -0.0382786504815158, -0.02346343035242171, -0.015845719190110685, -0.01141628299592412, -0.008614909256721148, -0.006731278790539363, -0.005404188868851634])
c_delta_low_a1 = np.array([-1.6948852145105775, -0.21221247998255421, -0.07616624070578837, -0.03925872051695478, -0.02392582928223419, -0.016100079392344924, -0.011571008417377016, -0.008715982312423876, -0.006800916617066832, -0.005454194160847692])


c_delta_fit_a01 =np.array([-1.568977772967628, -0.19475279777907417, -0.07363771574091516, -0.038435372425738024, -0.02355559963689302, -0.01590164802109939, -0.011452169474068796, -0.008639127827336779, -0.006748326723027276, -0.005416614749265136])

c_delta_low_a01 =np.array([-1.8238920248222712, -0.2031865221397311, -0.07555613929071114, -0.03915330080417334, -0.023898721428849967, -0.016091517409222433, -0.011568031823117053, -0.008714955311006634, -0.006800632309023058, -0.00545420343769365])


effective_fit_a1 = np.array([-1.5688546042838425, -0.19670531428346294, -0.07326681416088832, -0.03821402142420993, -0.023433193291566567, -0.015829150743229548, -0.011406224348320393, -0.008608345615357393, -0.006726759147568373, -0.00540094470125041])

effective_low_a1 =  np.array([-1.7896896096317505, -0.20940409076501965, -0.07599445461892174, -0.03922686855730717, -0.023916429290693486, -0.016096453236968955, -0.011569339676498203, -0.008715111562196398, -0.006800417486374499, -0.005453886569739552])

effective_fit_a01 =np.array([-1.5690337043452018, -0.19334354574311874, -0.07328466990657034, -0.038299369134620065, -0.023489748764404794, -0.01586494863659027, -0.011429676214902429, -0.008624363636045018, -0.0067381211010797415, -0.005409269124356797])
effective_low_a01 = np.array([-1.8216175551970082, -0.20310414829509682, -0.0755374617256166, -0.03914640228686039, -0.023895510003058007, -0.016089824475784553, -0.011567082492547343, -0.008714418163435766, -0.006800346545787761, -0.005454083930089837])
s=10

#plt.scatter(energy_levels,data,label='data')
plt.scatter(energy_levels,abs(data-coulomb)/data,label='Coulomb',s=s)
plt.scatter(energy_levels,abs(data-delta_1st_order_fit)/data,label='1st order $\delta$',s=s)
plt.scatter(energy_levels,abs(data-c_delta_fit_a1)/data,label='$c\delta a^2$, a=1',s=s)
plt.scatter(energy_levels,abs(data-c_delta_fit_a01)/data,label='$c\delta a^2$, a=0.1',s=s)
plt.scatter(energy_levels,abs(data-effective_fit_a1)/data,label='effective, a=1',s=s)
plt.scatter(energy_levels,abs(data-effective_fit_a01)/data,label='effective, a=0.1',s=s)
plt.legend()
plt.ylabel(r'$\Delta$E/E')
plt.xlabel('Energy levels')
plt.title("Energy with different model minimizing chi-square")

plt.savefig("pictures/final_results/graph_fit.png")
plt.close()

plt.scatter(energy_levels,abs(data-coulomb)/data,label=r'Coulomb',s=s)
plt.scatter(energy_levels,abs(data-delta_1st_order_low)/data,label=r'1st order $\delta$',s=s)
plt.scatter(energy_levels,abs(data-c_delta_low_a1)/data,label=r'$c\delta a^2$, a=1',s=s)
plt.scatter(energy_levels,abs(data-c_delta_low_a01)/data,label=r'$c\delta a^2$, a=0.1',s=s)
plt.scatter(energy_levels,abs(data-effective_low_a1)/data,label=r'effective, a=1',s=s)
plt.scatter(energy_levels,abs(data-effective_low_a01)/data,label=r'effective, a=0.1',s=s)
plt.legend()
plt.ylabel(r'$\Delta$E/E')
plt.xlabel('Energy levels')
plt.title("Energy with different model matching low energies with data")
plt.savefig("pictures/final_results/graph_lowenergy.png")

plt.close()