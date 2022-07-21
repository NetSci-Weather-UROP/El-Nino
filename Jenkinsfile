pipeline {
    agent any

	environment {
		WAYLAND_DISPLAY = 'wayland-0'
		DISPLAY = ':0'
	}

    stages {
		stage('Setup') {
			steps {
				step([$class: 'GitHubSetCommitStatusBuilder'])
				sh '[ -e air.sig995 ] || ln -s ~/air.sig995 .'
				sh '[ -e temp_data_1948_2021.npy ] || ln -s ~/temp_data_1948_2021.npy .'
				sh 'sed -i "s/plt.show()//g" nino.py'
			}
		}
		stage('Run') {
			parallel {
				stage('Run Python') {
					steps {
						sh 'python3 ./nino.py'
					}
				}
				stage('Run Julia') {
					steps {
						sh 'julia -t 1 -O 3 -C native ./main.jl'
					}
				}
			}
		}
    }
	post {
		 always {
		 	step([$class: 'GitHubCommitStatusSetter'])
		 }
	}
}
