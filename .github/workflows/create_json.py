from git import Repo
import sys, os, json
import argparse

def create_json (filename,
                 url_template,
                 include=[],
                 latest_version="dev") :
  """
  Create json file.
   
  ``url_template`` should be such as 
  ``https://mysite.org/``
  """
  repo = Repo( search_parent_directories=True )
  versions = [branch.name for branch in repo.branches]
  print ("Found corresponding branches on repo: {}".format (versions))
  json_list = []
  
  for version in versions:
    if version in include :
      if version==latest_version :
        version = "latest"
      version_dict = {
                      "version":version,
                      "url":os.path.join (url_template, version),
                      }
      json_list.append (version_dict)

  with open(filename, 'w') as fp:
    json.dump(json_list, fp)

if __name__ == "__main__" :

  CLI = argparse.ArgumentParser()
  CLI.add_argument(
    "--filename",  
    nargs="?", 
    type=str,
    default="switcher.json", 
    )
  CLI.add_argument(
    "--url",  
    nargs="?", 
    type=str,
    default="url_not_specified", 
    )
  CLI.add_argument(
    "--branches",  
    nargs="*", 
    type=str,
    default=["main", "dev"], 
    )
  CLI.add_argument(
    "--latest",  
    nargs="?", 
    type=str,
    default="dev", 
    )
  args = CLI.parse_args()
  print (vars(args).values())
  create_json (*vars(args).values())
