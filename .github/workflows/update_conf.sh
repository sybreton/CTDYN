# Update conf.py file in order to display version
# on the lower left 
cat >> docs/source/conf.py <<'EOF'

try:
   html_context
except NameError:
   html_context = dict()
html_context['display_lower_left'] = True

# set version
from git import Repo
repo = Repo( search_parent_directories=True )
if repo.active_branch.name=="dev" :
  current_version = "latest"
elif repo.active_branch.name=="main" :
  current_version = "stable"
else :
  current_version = repo.active_branch.name

html_context['current_version'] = current_version
html_context['version'] = current_version

versions = [branch.name for branch in repo.branches]
for version in versions:
   html_context['versions'].append( "_build/{}".format (version))

