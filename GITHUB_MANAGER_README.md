# GitHub Manager - Comprehensive GitHub Operations Tool

A standalone Python script that provides a complete toolkit for managing GitHub repositories, issues, pull requests, releases, and more through the GitHub API.

## Features

- **Repository Management**: Create, clone, fork, list repositories
- **Issue Management**: Create, list, close issues
- **Pull Request Management**: Create, list, merge pull requests  
- **Release Management**: Create and list releases
- **File Operations**: Upload, download, and update files
- **Workflow Management**: List and trigger workflows
- **User Management**: Get user information and search repositories
- **Statistics**: Get detailed repository statistics

## Installation

1. Install dependencies:
```bash
pip install -r requirements-github.txt
```

2. Set up your GitHub token:
```bash
export GITHUB_TOKEN="your_github_personal_access_token"
```

## Usage

### Get Help
```bash
python github_manager.py --help
```

### Repository Operations

**Create a repository:**
```bash
python github_manager.py repo create my-new-repo --description "My awesome project" --private
```

**List your repositories:**
```bash
python github_manager.py repo list
```

**List another user's repositories:**
```bash
python github_manager.py repo list --user username
```

**Clone a repository:**
```bash
python github_manager.py repo clone owner/repository-name --path ./local-folder
```

**Fork a repository:**
```bash
python github_manager.py repo fork owner/repository-name
```

**Get repository statistics:**
```bash
python github_manager.py repo stats owner/repository-name
```

### Issue Operations

**Create an issue:**
```bash
python github_manager.py issue create owner/repo "Bug report" --body "Description of the bug" --labels bug urgent
```

**List issues:**
```bash
python github_manager.py issue list owner/repo --state open
```

**Close an issue:**
```bash
python github_manager.py issue close owner/repo 123 --comment "Fixed in latest commit"
```

### Pull Request Operations

**Create a pull request:**
```bash
python github_manager.py pr create owner/repo "Add new feature" feature-branch --body "Description of changes" --base main
```

**List pull requests:**
```bash
python github_manager.py pr list owner/repo --state open
```

**Merge a pull request:**
```bash
python github_manager.py pr merge owner/repo 456 --message "Merge feature branch"
```

### Release Operations

**Create a release:**
```bash
python github_manager.py release create owner/repo v1.0.0 "Version 1.0.0" --body "Release notes"
```

**List releases:**
```bash
python github_manager.py release list owner/repo
```

### User and Search Operations

**Get user information:**
```bash
python github_manager.py user --username octocat
```

**Search repositories:**
```bash
python github_manager.py search "machine learning python" --sort stars --order desc
```

## Authentication

The script requires a GitHub Personal Access Token for authentication. You can:

1. **Set environment variable (recommended):**
```bash
export GITHUB_TOKEN="ghp_xxxxxxxxxxxxxxxxxxxx"
```

2. **Pass token as argument:**
```bash
python github_manager.py --token "ghp_xxxxxxxxxxxxxxxxxxxx" repo list
```

### Creating a GitHub Token

1. Go to GitHub Settings > Developer settings > Personal access tokens
2. Click "Generate new token"
3. Select appropriate scopes (repo, workflow, user, etc.)
4. Copy the generated token

## Programmatic Usage

You can also use the GitHubManager class directly in your Python scripts:

```python
from github_manager import GitHubManager

# Initialize manager
manager = GitHubManager("your_token_here")

# Create a repository
repo = manager.create_repository("my-repo", "A test repository")

# Create an issue
issue = manager.create_issue("owner/repo", "Bug report", "Description")

# List repositories
repos = manager.list_repositories()
```

## Error Handling

The script includes comprehensive error handling and logging. All operations will log their status and any errors encountered.

## Requirements

- Python 3.6+
- PyGithub library
- Valid GitHub Personal Access Token

## Security Notes

- Never commit your GitHub token to version control
- Use environment variables for token storage
- Limit token permissions to only what you need
- Regularly rotate your tokens for security