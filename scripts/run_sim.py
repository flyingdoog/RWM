import os
import time
device_number = 7
windows = True

default_settings={
	'dataset': None,
	'alpha' :0.9,
	'lamda' : 0.7,
	'theta' : 0.9,
	'epsilon_matrix':0.01,
	'step' : 1,
	'min_cmty_size' : 5,
	'max_cmty_size' :200,
	'ite1st':-1,
	'mute':1
}

datasets = ['rm'] #RM,dblp,brainNet

steps = {
	'rm':1,
	'brainnet':1,
}

default_alphas = {
	'rm':0.7,
	'brainnet':0.9,
}

default_lambdas = {
	'rm':0.9,
	'brainnet':0.7,
}

cmty_sizes = {
	'rm':(5,200),
	'brainnet':(5,200),
}

alphas = [(i)/10.0 for i in range(1,10)]
thetas = [(i)/10.0 for i in range(1,10)]
lamdas = [(i)/10.0 for i in range(1,10)]
ite1sts = [i for i in range(0,10)]
print(thetas)

# 
try:
	if windows:
		os.system('del /Q flags\\*')
	else:
		os.system('rm -rf flags')
		os.system('mkdir flags')

except:
	pass

def get_available_device():

	while True:
		for i in range(device_number):
			if not str(i) in os.listdir('flags'):
				if windows:
					os.system("echo nul > flags\\"+str(i))
				else:
					os.system('touch flags/'+str(i))
				time.sleep(1)
				return i
		print('waiting')
		time.sleep(120)

def run(settings,device_id):
	if windows:
		cmd = 'rwm.exe '
	else:
		cmd = '(./RWM '
	cmd += ' -d '+settings['dataset']+' '
	cmd += ' -a '+str(settings['alpha'])+' '
	cmd += ' -l '+str(settings['lamda'])+' '
	cmd += ' -t '+str(settings['theta'])+' '
	cmd += ' -e '+str(settings['epsilon_matrix'])+' '
	cmd += ' -s '+str(settings['step'])+' '
	cmd += ' -min '+str(settings['min_cmty_size'])+' '
	cmd += ' -max '+str(settings['max_cmty_size'])+ ' '
	cmd += ' -ite1st '+str(settings['ite1st'])
	cmd += ' -m '+str(settings['mute'])
	if windows:
		cmd +=  '&& del /Q flags\\'+str(device_id)+''
	else:
		cmd += '; rm flags/'+str(device_id)+')&'
	print(cmd)
	os.system(cmd)

def run_whole(settings):
	device_id = get_available_device()
	run(ts,device_id)

def run_seq(settings):
	if windows:
		cmd = 'rwm.exe '
	else:
		cmd = './RWM '
	cmd += ' -d '+settings['dataset']+' '
	cmd += ' -a '+str(settings['alpha'])+' '
	cmd += ' -l '+str(settings['lamda'])+' '
	cmd += ' -t '+str(settings['theta'])+' '
	cmd += ' -e '+str(settings['epsilon_matrix'])+' '
	cmd += ' -s '+str(settings['step'])+' '
	cmd += ' -min '+str(settings['min_cmty_size'])+' '
	cmd += ' -max '+str(settings['max_cmty_size'])+ ' '
	cmd += ' -ite1st '+str(settings['ite1st'])
	cmd += ' -m '+str(settings['mute'])
	print(cmd)
	os.system(cmd)


for dataset in datasets:
	settings = default_settings.copy()
	settings['dataset'] = dataset
	settings['min_cmty_size'] = cmty_sizes[dataset][0]
	settings['max_cmty_size'] = cmty_sizes[dataset][1]
	settings['step'] = steps[dataset]
	settings['alpha'] = default_alphas[dataset]
	settings['lamda'] = default_lambdas[dataset]


	for alpha in alphas:
		ts = settings.copy();
		ts['alpha']=alpha
		run_seq(ts)

	for lamda in lamdas:
		ts = settings.copy();
		ts['lamda']=lamda
		run_seq(ts)
	
	for theta in thetas:
		ts = settings.copy();
		ts['theta']=theta
		run_seq(ts)


	for ite1st in ite1sts:
		ts = settings.copy();
		ts['ite1st']=ite1st
		run_seq(ts)

