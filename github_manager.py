#!/usr/bin/env python3
"""
GitHub Manager - A comprehensive standalone Python script for GitHub operations
==============================================================================

This script provides a complete toolkit for managing GitHub repositories,
issues, pull requests, releases, and more through the GitHub API.

Features:
- Repository management (create, clone, fork, etc.)
- Issue management (create, update, close, list)
- Pull request management (create, merge, review)
- Release management (create, update, list)
- File operations (upload, download, edit)
- Workflow management
- Organization management
- User management

Requirements:
- PyGithub library (pip install PyGithub)
- GitHub personal access token

Usage:
    python github_manager.py --help
"""

import argparse
import json
import os
import sys
from datetime import datetime
from typing import Dict, List, Optional, Any
import subprocess
import logging

try:
    from github import Github, GithubException
    from github.Repository import Repository
    from github.Issue import Issue
    from github.PullRequest import PullRequest
    from github.GitRelease import GitRelease as Release
except ImportError:
    print("PyGithub library not found. Please install it with: pip install PyGithub")
    sys.exit(1)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GitHubManager:
    """Comprehensive GitHub management class"""
    
    def __init__(self, token: str = None):
        """Initialize GitHub manager with authentication token"""
        self.token = token or os.environ.get('GITHUB_TOKEN')
        if not self.token:
            raise ValueError("GitHub token is required. Set GITHUB_TOKEN environment variable or pass token parameter.")
        
        try:
            self.github = Github(self.token)
            self.user = self.github.get_user()
            logger.info(f"Successfully authenticated as: {self.user.login}")
        except GithubException as e:
            logger.error(f"Authentication failed: {e}")
            raise
    
    def get_repository(self, repo_name: str) -> Repository:
        """Get repository object"""
        try:
            return self.github.get_repo(repo_name)
        except GithubException as e:
            logger.error(f"Failed to get repository {repo_name}: {e}")
            raise
    
    # Repository Operations
    def create_repository(self, name: str, description: str = "", private: bool = False, 
                         has_issues: bool = True, has_projects: bool = True, 
                         has_wiki: bool = True) -> Repository:
        """Create a new repository"""
        try:
            repo = self.user.create_repo(
                name=name,
                description=description,
                private=private,
                has_issues=has_issues,
                has_projects=has_projects,
                has_wiki=has_wiki
            )
            logger.info(f"Successfully created repository: {repo.full_name}")
            return repo
        except GithubException as e:
            logger.error(f"Failed to create repository: {e}")
            raise
    
    def clone_repository(self, repo_name: str, local_path: str = None) -> str:
        """Clone a repository to local machine"""
        try:
            repo = self.get_repository(repo_name)
            clone_url = repo.clone_url
            
            if local_path is None:
                local_path = repo.name
            
            cmd = ["git", "clone", clone_url, local_path]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info(f"Successfully cloned {repo_name} to {local_path}")
                return local_path
            else:
                logger.error(f"Failed to clone repository: {result.stderr}")
                raise Exception(result.stderr)
        except Exception as e:
            logger.error(f"Error cloning repository: {e}")
            raise
    
    def fork_repository(self, repo_name: str) -> Repository:
        """Fork a repository"""
        try:
            repo = self.get_repository(repo_name)
            forked_repo = self.user.create_fork(repo)
            logger.info(f"Successfully forked {repo_name} to {forked_repo.full_name}")
            return forked_repo
        except GithubException as e:
            logger.error(f"Failed to fork repository: {e}")
            raise
    
    def list_repositories(self, user: str = None) -> List[Repository]:
        """List repositories for a user (default: authenticated user)"""
        try:
            if user:
                github_user = self.github.get_user(user)
                repos = github_user.get_repos()
            else:
                repos = self.user.get_repos()
            
            repo_list = list(repos)
            logger.info(f"Found {len(repo_list)} repositories")
            return repo_list
        except GithubException as e:
            logger.error(f"Failed to list repositories: {e}")
            raise
    
    # Issue Operations
    def create_issue(self, repo_name: str, title: str, body: str = "", 
                    labels: List[str] = None, assignees: List[str] = None) -> Issue:
        """Create a new issue"""
        try:
            repo = self.get_repository(repo_name)
            issue = repo.create_issue(
                title=title,
                body=body,
                labels=labels or [],
                assignees=assignees or []
            )
            logger.info(f"Successfully created issue #{issue.number}: {title}")
            return issue
        except GithubException as e:
            logger.error(f"Failed to create issue: {e}")
            raise
    
    def list_issues(self, repo_name: str, state: str = "open") -> List[Issue]:
        """List issues in a repository"""
        try:
            repo = self.get_repository(repo_name)
            issues = repo.get_issues(state=state)
            issue_list = list(issues)
            logger.info(f"Found {len(issue_list)} {state} issues")
            return issue_list
        except GithubException as e:
            logger.error(f"Failed to list issues: {e}")
            raise
    
    def close_issue(self, repo_name: str, issue_number: int, comment: str = "") -> Issue:
        """Close an issue"""
        try:
            repo = self.get_repository(repo_name)
            issue = repo.get_issue(issue_number)
            
            if comment:
                issue.create_comment(comment)
            
            issue.edit(state="closed")
            logger.info(f"Successfully closed issue #{issue_number}")
            return issue
        except GithubException as e:
            logger.error(f"Failed to close issue: {e}")
            raise
    
    # Pull Request Operations
    def create_pull_request(self, repo_name: str, title: str, body: str, 
                          head: str, base: str = "main") -> PullRequest:
        """Create a new pull request"""
        try:
            repo = self.get_repository(repo_name)
            pr = repo.create_pull(
                title=title,
                body=body,
                head=head,
                base=base
            )
            logger.info(f"Successfully created pull request #{pr.number}: {title}")
            return pr
        except GithubException as e:
            logger.error(f"Failed to create pull request: {e}")
            raise
    
    def list_pull_requests(self, repo_name: str, state: str = "open") -> List[PullRequest]:
        """List pull requests in a repository"""
        try:
            repo = self.get_repository(repo_name)
            prs = repo.get_pulls(state=state)
            pr_list = list(prs)
            logger.info(f"Found {len(pr_list)} {state} pull requests")
            return pr_list
        except GithubException as e:
            logger.error(f"Failed to list pull requests: {e}")
            raise
    
    def merge_pull_request(self, repo_name: str, pr_number: int, 
                          commit_message: str = "") -> bool:
        """Merge a pull request"""
        try:
            repo = self.get_repository(repo_name)
            pr = repo.get_pull(pr_number)
            
            merge_result = pr.merge(commit_message=commit_message)
            logger.info(f"Successfully merged pull request #{pr_number}")
            return merge_result.merged
        except GithubException as e:
            logger.error(f"Failed to merge pull request: {e}")
            raise
    
    # Release Operations
    def create_release(self, repo_name: str, tag_name: str, name: str, 
                      body: str = "", draft: bool = False, prerelease: bool = False) -> Release:
        """Create a new release"""
        try:
            repo = self.get_repository(repo_name)
            release = repo.create_git_release(
                tag=tag_name,
                name=name,
                message=body,
                draft=draft,
                prerelease=prerelease
            )
            logger.info(f"Successfully created release: {name}")
            return release
        except GithubException as e:
            logger.error(f"Failed to create release: {e}")
            raise
    
    def list_releases(self, repo_name: str) -> List[Release]:
        """List releases in a repository"""
        try:
            repo = self.get_repository(repo_name)
            releases = repo.get_releases()
            release_list = list(releases)
            logger.info(f"Found {len(release_list)} releases")
            return release_list
        except GithubException as e:
            logger.error(f"Failed to list releases: {e}")
            raise
    
    # File Operations
    def upload_file(self, repo_name: str, file_path: str, content: str, 
                   commit_message: str, branch: str = "main") -> bool:
        """Upload a file to repository"""
        try:
            repo = self.get_repository(repo_name)
            repo.create_file(
                path=file_path,
                message=commit_message,
                content=content,
                branch=branch
            )
            logger.info(f"Successfully uploaded file: {file_path}")
            return True
        except GithubException as e:
            logger.error(f"Failed to upload file: {e}")
            raise
    
    def download_file(self, repo_name: str, file_path: str, branch: str = "main") -> str:
        """Download a file from repository"""
        try:
            repo = self.get_repository(repo_name)
            file_content = repo.get_contents(file_path, ref=branch)
            content = file_content.decoded_content.decode('utf-8')
            logger.info(f"Successfully downloaded file: {file_path}")
            return content
        except GithubException as e:
            logger.error(f"Failed to download file: {e}")
            raise
    
    def update_file(self, repo_name: str, file_path: str, content: str, 
                   commit_message: str, branch: str = "main") -> bool:
        """Update an existing file in repository"""
        try:
            repo = self.get_repository(repo_name)
            file_content = repo.get_contents(file_path, ref=branch)
            
            repo.update_file(
                path=file_path,
                message=commit_message,
                content=content,
                sha=file_content.sha,
                branch=branch
            )
            logger.info(f"Successfully updated file: {file_path}")
            return True
        except GithubException as e:
            logger.error(f"Failed to update file: {e}")
            raise
    
    # Workflow Operations
    def list_workflows(self, repo_name: str) -> List[Dict]:
        """List workflows in a repository"""
        try:
            repo = self.get_repository(repo_name)
            workflows = repo.get_workflows()
            workflow_list = []
            
            for workflow in workflows:
                workflow_list.append({
                    'id': workflow.id,
                    'name': workflow.name,
                    'path': workflow.path,
                    'state': workflow.state,
                    'created_at': workflow.created_at,
                    'updated_at': workflow.updated_at
                })
            
            logger.info(f"Found {len(workflow_list)} workflows")
            return workflow_list
        except GithubException as e:
            logger.error(f"Failed to list workflows: {e}")
            raise
    
    def trigger_workflow(self, repo_name: str, workflow_id: str, 
                        ref: str = "main", inputs: Dict = None) -> bool:
        """Trigger a workflow run"""
        try:
            repo = self.get_repository(repo_name)
            workflow = repo.get_workflow(workflow_id)
            
            success = workflow.create_dispatch(ref=ref, inputs=inputs or {})
            logger.info(f"Successfully triggered workflow: {workflow.name}")
            return success
        except GithubException as e:
            logger.error(f"Failed to trigger workflow: {e}")
            raise
    
    # Utility Methods
    def get_user_info(self, username: str = None) -> Dict:
        """Get user information"""
        try:
            if username:
                user = self.github.get_user(username)
            else:
                user = self.user
            
            return {
                'login': user.login,
                'name': user.name,
                'email': user.email,
                'company': user.company,
                'location': user.location,
                'bio': user.bio,
                'public_repos': user.public_repos,
                'followers': user.followers,
                'following': user.following,
                'created_at': user.created_at,
                'updated_at': user.updated_at
            }
        except GithubException as e:
            logger.error(f"Failed to get user info: {e}")
            raise
    
    def search_repositories(self, query: str, sort: str = "stars", 
                          order: str = "desc") -> List[Repository]:
        """Search repositories"""
        try:
            repos = self.github.search_repositories(query, sort=sort, order=order)
            repo_list = list(repos)
            logger.info(f"Found {len(repo_list)} repositories matching '{query}'")
            return repo_list
        except GithubException as e:
            logger.error(f"Failed to search repositories: {e}")
            raise
    
    def get_repository_stats(self, repo_name: str) -> Dict:
        """Get repository statistics"""
        try:
            repo = self.get_repository(repo_name)
            
            return {
                'name': repo.name,
                'full_name': repo.full_name,
                'description': repo.description,
                'language': repo.language,
                'stars': repo.stargazers_count,
                'forks': repo.forks_count,
                'watchers': repo.watchers_count,
                'issues': repo.open_issues_count,
                'size': repo.size,
                'created_at': repo.created_at,
                'updated_at': repo.updated_at,
                'pushed_at': repo.pushed_at,
                'default_branch': repo.default_branch,
                'license': repo.license.name if repo.license else None,
                'private': repo.private,
                'archived': repo.archived,
                'disabled': repo.disabled
            }
        except GithubException as e:
            logger.error(f"Failed to get repository stats: {e}")
            raise

def main():
    """Main function with CLI interface"""
    parser = argparse.ArgumentParser(description='GitHub Manager - Comprehensive GitHub operations tool')
    parser.add_argument('--token', help='GitHub personal access token')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Repository commands
    repo_parser = subparsers.add_parser('repo', help='Repository operations')
    repo_subparsers = repo_parser.add_subparsers(dest='repo_action')
    
    # Create repository
    create_repo_parser = repo_subparsers.add_parser('create', help='Create a new repository')
    create_repo_parser.add_argument('name', help='Repository name')
    create_repo_parser.add_argument('--description', help='Repository description')
    create_repo_parser.add_argument('--private', action='store_true', help='Make repository private')
    
    # List repositories
    list_repos_parser = repo_subparsers.add_parser('list', help='List repositories')
    list_repos_parser.add_argument('--user', help='Username (default: authenticated user)')
    
    # Clone repository
    clone_parser = repo_subparsers.add_parser('clone', help='Clone a repository')
    clone_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    clone_parser.add_argument('--path', help='Local path to clone to')
    
    # Fork repository
    fork_parser = repo_subparsers.add_parser('fork', help='Fork a repository')
    fork_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    
    # Repository stats
    stats_parser = repo_subparsers.add_parser('stats', help='Get repository statistics')
    stats_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    
    # Issue commands
    issue_parser = subparsers.add_parser('issue', help='Issue operations')
    issue_subparsers = issue_parser.add_subparsers(dest='issue_action')
    
    # Create issue
    create_issue_parser = issue_subparsers.add_parser('create', help='Create a new issue')
    create_issue_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    create_issue_parser.add_argument('title', help='Issue title')
    create_issue_parser.add_argument('--body', help='Issue body')
    create_issue_parser.add_argument('--labels', nargs='+', help='Issue labels')
    
    # List issues
    list_issues_parser = issue_subparsers.add_parser('list', help='List issues')
    list_issues_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    list_issues_parser.add_argument('--state', choices=['open', 'closed', 'all'], default='open', help='Issue state')
    
    # Close issue
    close_issue_parser = issue_subparsers.add_parser('close', help='Close an issue')
    close_issue_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    close_issue_parser.add_argument('issue_number', type=int, help='Issue number')
    close_issue_parser.add_argument('--comment', help='Closing comment')
    
    # Pull request commands
    pr_parser = subparsers.add_parser('pr', help='Pull request operations')
    pr_subparsers = pr_parser.add_subparsers(dest='pr_action')
    
    # Create pull request
    create_pr_parser = pr_subparsers.add_parser('create', help='Create a new pull request')
    create_pr_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    create_pr_parser.add_argument('title', help='Pull request title')
    create_pr_parser.add_argument('head', help='Head branch')
    create_pr_parser.add_argument('--body', help='Pull request body')
    create_pr_parser.add_argument('--base', default='main', help='Base branch (default: main)')
    
    # List pull requests
    list_prs_parser = pr_subparsers.add_parser('list', help='List pull requests')
    list_prs_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    list_prs_parser.add_argument('--state', choices=['open', 'closed', 'all'], default='open', help='PR state')
    
    # Merge pull request
    merge_pr_parser = pr_subparsers.add_parser('merge', help='Merge a pull request')
    merge_pr_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    merge_pr_parser.add_argument('pr_number', type=int, help='Pull request number')
    merge_pr_parser.add_argument('--message', help='Merge commit message')
    
    # Release commands
    release_parser = subparsers.add_parser('release', help='Release operations')
    release_subparsers = release_parser.add_subparsers(dest='release_action')
    
    # Create release
    create_release_parser = release_subparsers.add_parser('create', help='Create a new release')
    create_release_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    create_release_parser.add_argument('tag_name', help='Release tag name')
    create_release_parser.add_argument('name', help='Release name')
    create_release_parser.add_argument('--body', help='Release body')
    create_release_parser.add_argument('--draft', action='store_true', help='Create as draft')
    create_release_parser.add_argument('--prerelease', action='store_true', help='Mark as prerelease')
    
    # List releases
    list_releases_parser = release_subparsers.add_parser('list', help='List releases')
    list_releases_parser.add_argument('repo_name', help='Repository name (owner/repo)')
    
    # User info
    user_parser = subparsers.add_parser('user', help='Get user information')
    user_parser.add_argument('--username', help='Username (default: authenticated user)')
    
    # Search repositories
    search_parser = subparsers.add_parser('search', help='Search repositories')
    search_parser.add_argument('query', help='Search query')
    search_parser.add_argument('--sort', choices=['stars', 'forks', 'updated'], default='stars', help='Sort by')
    search_parser.add_argument('--order', choices=['asc', 'desc'], default='desc', help='Sort order')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    try:
        # Initialize GitHub manager
        github_manager = GitHubManager(args.token)
        
        # Execute commands
        if args.command == 'repo':
            if args.repo_action == 'create':
                repo = github_manager.create_repository(
                    args.name, 
                    args.description or "", 
                    args.private
                )
                print(f"Created repository: {repo.html_url}")
            
            elif args.repo_action == 'list':
                repos = github_manager.list_repositories(args.user)
                for repo in repos:
                    print(f"{repo.full_name} - {repo.description}")
            
            elif args.repo_action == 'clone':
                path = github_manager.clone_repository(args.repo_name, args.path)
                print(f"Cloned repository to: {path}")
            
            elif args.repo_action == 'fork':
                repo = github_manager.fork_repository(args.repo_name)
                print(f"Forked repository: {repo.html_url}")
            
            elif args.repo_action == 'stats':
                stats = github_manager.get_repository_stats(args.repo_name)
                print(json.dumps(stats, indent=2, default=str))
        
        elif args.command == 'issue':
            if args.issue_action == 'create':
                issue = github_manager.create_issue(
                    args.repo_name, 
                    args.title, 
                    args.body or "", 
                    args.labels
                )
                print(f"Created issue: {issue.html_url}")
            
            elif args.issue_action == 'list':
                issues = github_manager.list_issues(args.repo_name, args.state)
                for issue in issues:
                    print(f"#{issue.number}: {issue.title}")
            
            elif args.issue_action == 'close':
                issue = github_manager.close_issue(args.repo_name, args.issue_number, args.comment)
                print(f"Closed issue #{issue.number}")
        
        elif args.command == 'pr':
            if args.pr_action == 'create':
                pr = github_manager.create_pull_request(
                    args.repo_name, 
                    args.title, 
                    args.body or "", 
                    args.head, 
                    args.base
                )
                print(f"Created pull request: {pr.html_url}")
            
            elif args.pr_action == 'list':
                prs = github_manager.list_pull_requests(args.repo_name, args.state)
                for pr in prs:
                    print(f"#{pr.number}: {pr.title}")
            
            elif args.pr_action == 'merge':
                success = github_manager.merge_pull_request(
                    args.repo_name, 
                    args.pr_number, 
                    args.message
                )
                print(f"Pull request merged: {success}")
        
        elif args.command == 'release':
            if args.release_action == 'create':
                release = github_manager.create_release(
                    args.repo_name, 
                    args.tag_name, 
                    args.name, 
                    args.body or "", 
                    args.draft, 
                    args.prerelease
                )
                print(f"Created release: {release.html_url}")
            
            elif args.release_action == 'list':
                releases = github_manager.list_releases(args.repo_name)
                for release in releases:
                    print(f"{release.tag_name}: {release.name}")
        
        elif args.command == 'user':
            user_info = github_manager.get_user_info(args.username)
            print(json.dumps(user_info, indent=2, default=str))
        
        elif args.command == 'search':
            repos = github_manager.search_repositories(args.query, args.sort, args.order)
            for repo in repos[:10]:  # Limit to first 10 results
                print(f"{repo.full_name} - ⭐ {repo.stargazers_count} - {repo.description}")
    
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()