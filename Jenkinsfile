pipeline {
    agent any
    stages {
        stage('Run') {
            steps {
				step([$class: 'GitHubSetCommitStatusBuilder'])
                sh 'exec ./nino.jl'
            }
        }
    }
	post {
		 always {
		 	step([$class: 'GitHubCommitStatusSetter'])
		 }
	}
}
