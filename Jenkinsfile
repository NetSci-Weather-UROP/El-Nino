pipeline {
    agent any
    stages {
		stage('Setup') {
			steps {
				step([$class: 'GitHubSetCommitStatusBuilder'])
				sh 'ln -s ~/air.sig995 .'
			}
		}
        stage('Run') {
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
