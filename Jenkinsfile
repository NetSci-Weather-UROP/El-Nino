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
				sh 'sed -i "s/plt.show()/pass/g" nino.py'
				sh 'rm CNW-plots/*'
			}
		}
		stage('Run') {
			parallel {
				stage('Run Python') {
					steps {
						sh 'python3 ./nino.py'
						archiveArtifacts artifacts: 'CNW-plots/*'
						sh 'sleep 2'
						sh 'echo >> wait'
					}
				}
				stage('Run Julia') {
					steps {
						sh 'mkfifo wait'
						sh 'cat wait'
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
