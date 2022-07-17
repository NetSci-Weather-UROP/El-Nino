pipeline {
    agent any
    stages {
		stage('Setup') {
			steps {
				step([$class: 'GitHubSetCommitStatusBuilder'])
				sh '[ -e air.sig995 ] || ln -s ~/air.sig995 .'
			}
		}
		stage('Run Python') {
			steps {
				sh 'python3 ./nino.jl'
			}
		}
        stage('Run Julia') {
            steps {
                sh 'julia -t 1 -O 3 -C native ./main.jl'
            }
        }
    }
	post {
		 always {
		 	step([$class: 'GitHubCommitStatusSetter'])
		 }
	}
}
