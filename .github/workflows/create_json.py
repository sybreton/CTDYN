from git import Repo
import sys, os, json

def create_json (filename,
                 url_template,
                 include=[]) :
  """
  Create json file.
   
  ``url_template`` should be such as 
  ``https://mysite.org/``
  """
  repo = Repo( search_parent_directories=True )
  
  versions = [branch.name for branch in repo.branches]
  
  json_list = []
  
  for version in versions:
    if version in include :
      version_dict = {
                      "version":version,
                      "url":os.path.join (url_template, version),
                      }
      json_list.append (version_dict)

  with open(filename, 'w') as fp:
    json.dump(json_list, fp)

if __name__ == "__main__" :
  create_json (sys.argv[0], sys.argv[1], sys.argv[2])
