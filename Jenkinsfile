pipeline {
    agent any
    stages {
        stage('Run') {
            steps {
				step([$class: 'GitHubSetCommitStatusBuilder'])
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
