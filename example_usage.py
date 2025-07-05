#!/usr/bin/env python3
"""
Example usage of the GitHubManager class
=========================================

This script demonstrates how to use the GitHubManager class
programmatically for various GitHub operations.

Before running this script:
1. Set your GitHub token: export GITHUB_TOKEN="your_token_here"
2. Install dependencies: pip3 install -r requirements-github.txt
"""

from github_manager import GitHubManager
import os

def example_usage():
    """Demonstrate various GitHub operations"""
    
    # Initialize the GitHub manager
    # Token can be passed directly or via environment variable
    github_token = os.environ.get('GITHUB_TOKEN')
    if not github_token:
        print("Please set your GITHUB_TOKEN environment variable")
        print("export GITHUB_TOKEN='your_github_token_here'")
        return
    
    try:
        # Initialize manager
        manager = GitHubManager(github_token)
        print(f"✅ Successfully authenticated as: {manager.user.login}")
        
        # Example 1: Get user information
        print("\n" + "="*50)
        print("1. Getting user information")
        print("="*50)
        user_info = manager.get_user_info()
        print(f"Name: {user_info['name']}")
        print(f"Username: {user_info['login']}")
        print(f"Public repos: {user_info['public_repos']}")
        print(f"Followers: {user_info['followers']}")
        
        # Example 2: List your repositories
        print("\n" + "="*50)
        print("2. Listing your repositories (first 5)")
        print("="*50)
        repos = manager.list_repositories()
        for i, repo in enumerate(repos[:5]):
            print(f"{i+1}. {repo.full_name} - {repo.description or 'No description'}")
        
        # Example 3: Search repositories
        print("\n" + "="*50)
        print("3. Searching for Python repositories (first 3)")
        print("="*50)
        search_results = manager.search_repositories("python machine learning", sort="stars")
        for i, repo in enumerate(search_results[:3]):
            print(f"{i+1}. {repo.full_name} - ⭐ {repo.stargazers_count}")
            print(f"   {repo.description or 'No description'}")
        
        # Example 4: Get repository statistics (for a public repo)
        print("\n" + "="*50)
        print("4. Getting repository statistics for 'octocat/Hello-World'")
        print("="*50)
        try:
            stats = manager.get_repository_stats("octocat/Hello-World")
            print(f"Stars: {stats['stars']}")
            print(f"Forks: {stats['forks']}")
            print(f"Language: {stats['language']}")
            print(f"Created: {stats['created_at']}")
        except Exception as e:
            print(f"Could not get stats: {e}")
        
        # Example 5: List issues for a repository (if you have access)
        print("\n" + "="*50)
        print("5. Example of listing issues (commented out - needs repo access)")
        print("="*50)
        print("# To list issues from a repository you have access to:")
        print("# issues = manager.list_issues('your-username/your-repo')")
        print("# for issue in issues[:3]:")
        print("#     print(f'#{issue.number}: {issue.title}')")
        
        print("\n✅ All examples completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    example_usage()