import gitlab
import os
import re
import operator
import uuid


def main():
    server_url="https://git.rwth-aachen.de/"
    project_name_with_namespace = "aia/MAIA/Solver"

    if os.environ.get('PRIVATE_TOKEN_FILE') is not None:
        with open(os.environ['PRIVATE_TOKEN_FILE']) as f:
            token = f.readlines()[0].rstrip()
    else:
        token = os.environ['PRIVATE_TOKEN']

    gl = gitlab.Gitlab(server_url, private_token=token)
    # projects = gl.projects.list()
    project = gl.projects.get(project_name_with_namespace)
    # print(project)

    # Get open project issues
    issues = project.issues.list(state='opened', all=True)

    general_todo_issue = ""
    todo_issues = []
    # Search for general TODOs issue
    for issue in issues:
        # print(issue)
        # print(issue.title)
        # print(issue.iid)
        title = issue.title
        if title == "TODOs":
            general_todo_issue = issue
            print("Found general TODOs issue: '{}'; iid: {}".format(title, issue.iid))

        if title.startswith("TODOs - "):
            issue_label = title.replace("TODOs - ", "")
            todo_issues.append([issue_label, issue])
            print("Found TODOs issue: '{}'; label: {}; iid: {}".format(title, issue_label, issue.iid))

    # if general_todo_issue != "":
    #     print("Issue already existing", general_todo_issue.iid)

    # Not required for pipeline! TODO check if called from CI pipeline or manual!
    # path_ext = str(uuid.uuid4())
    # git_path="/tmp/MAIA-Solver_"+path_ext
    # repo_url="git@git.rwth-aachen.de:aia/MAIA/Solver.git"

    # os.system("git clone " + repo_url + " " + git_path)

    git_path = os.getcwd()
    # grep for keywords: TODO, FIXME, BUG, HACK, REVIEW
    todo_grep=os.popen('grep -iRn "\<TODO\|\<FIXME\|\<BUG\|\<HACK\|\<REVIEW" --include \*.cpp --include \*.h '+git_path+'/src/*').read()

    general_todos_description = ""
    todos_list = []
    # labels = set({""})
    labels = set()
    # print(type(labels))

    for line in todo_grep.splitlines():
        line = line.replace(git_path, "")
        line_split = re.findall(".*?:", line)
        filename = line_split[0].replace(":", "")
        linenumber = line_split[1].replace(":", "")

        line_code = re.findall(linenumber+":.*", line)[0].replace(linenumber+":", "")
        # print(line_code)
        code_labels = re.findall("labels:[a-zA-Z,]*", line_code)
        if len(code_labels) > 0:
            code_labels = code_labels[0].replace("labels:", "").split(",")
            # print(code_labels)
            # for l in code_labels:
            #     print("add", l)
            #     labels.add(l)
            #     print(labels)

            # labels.update(x for x in code_labels)
            labels.update(code_labels)

        first_line = int(linenumber)-2
        last_line = int(linenumber)+2
        sed_command="sed -n '"+str(first_line)+","+str(last_line)+"p;"+str(last_line+1)+"q' "+git_path+filename
        code_snippet = os.popen(sed_command).read()

        todos_list.append([filename, linenumber, code_snippet, code_labels, 0])

    # only update existing issues! -> todos with labels that does not match go to general todos issue
    print("Labels found in code: {}".format(labels))
    labels_tmp = labels
    labels = set()
    for label in labels_tmp:
        issue = [x for x in todo_issues if x[0]==label]
        if len(issue) == 1:
            labels.add(label)
            # issue = issue[0]
            print("Issue exists", issue)
        elif len(issue) == 0:
            print("Issue for label '{}' does not exist.".format(label))
        else:
            print("ERROR: multiple issues found for label {}".format(label))
            exit(1)
    print("Labels with a corresponding issue: {}".format(labels))

    # Sort by filename and by linenumber
    sorted(todos_list, key=operator.itemgetter(1,2), reverse=False)

    # TODO check for issues that exist in gitlab but have no corresponding label left in the code!
    for label in labels:
        todos_description = "# TODO/FIXME/BUG/HACK/REVIEW keywords found in source code with label {}\n".format(label)
        todos_description += "[[_TOC_]]\n\n"

        current_file = ""
        # Assemble issue description
        for todo in todos_list:
            # check if label is present
            if label not in todo[3]:
                continue
            # increase count of how many times this todo is part of an issue
            todo[4] += 1

            if todo[0] != current_file:
                # todo_count_file = sum(t[0].count(todo[0]) for t in todos_list)
                todos_description += "## {}\n".format(todo[0])
                # + " - " + str(todo_count_file) + "\n"
                current_file = todo[0]

            location = "{}#L{}".format(todo[0], todo[1])
            todos_description += "[{}]({}) - labels: {}\n".format(location, location, ', '.join(todo[3]))
            todos_description += "``` c++\n"
            todos_description += todo[2]
            todos_description += "```\n"
        # print(todos_description)

        issue = [x for x in todo_issues if x[0]==label]
        if len(issue) == 1:
            issue = issue[0]
            print("Checking issue {}".format(issue))
            descr_changed = check_description_change(issue[1].description, todos_description)
            if descr_changed:
                print("Issue description changed.", descr_changed)
                print("Updating issue", issue[1].iid)
                gl_issue = project.issues.get(issue[1].iid, lazy=True)
                gl_issue.description = todos_description
                gl_issue.save()
            else:
                print("Issue description for issue {} with label '{}' has not changed.".format(issue[1].iid, label))
        else:
            print("ERROR: label {}".format(label))
            exit(1)


    todo_count_total = len([t for t in todos_list if t[4]==0])
    general_todos_description += "# TODO/FIXME/BUG/HACK/REVIEW keywords found in source code - $`\sum`$" + str(todo_count_total) + "\n"
    # table of contents
    general_todos_description += "[[_TOC_]]\n\n"

    current_file = ""
    # Assemble issue description
    for todo in todos_list:
        # Skip todos already included in a labeled issue
        if todo[4] > 0:
            # print("Skipping TODO {}".format(todo[4]))
            continue

        if todo[0] != current_file:
            todo_count_file = sum(t[0].count(todo[0]) for t in todos_list if t[4]==0)
            general_todos_description += "## " + todo[0] + " - " + str(todo_count_file) + "\n"
            current_file = todo[0]

        location = todo[0]+"#L"+todo[1]
        general_todos_description += "["+location+"]("+location+")\n"
        general_todos_description += "``` c++\n"
        general_todos_description += todo[2]
        general_todos_description += "```\n"

    # print(general_todos_description)

    # issue_details = {
    # 'title': "TODOs",
    # 'description': general_todos_description
    # }

    if general_todo_issue == "":
        print("ERROR: General TODO issue should already exist.")
        exit(1)
        # print("Creating new issue")
        # new_issue = project.issues.create(issue_details)
        # print("New issue id {}".format(new_issue.iid))
    else:
        descr_changed = check_description_change(general_todo_issue.description, general_todos_description)
        if descr_changed:
            print("Update issue", general_todo_issue.iid)
            gl_issue = project.issues.get(general_todo_issue.iid, lazy=True)
            gl_issue.description = general_todos_description
            gl_issue.save()
        else:
            print("Issue description for general TODOs issue {} not changed".format(general_todo_issue.iid))

def check_description_change(old_descr, new_descr):
    newfile = '/tmp/new_description_{}.txt'.format(uuid.uuid4())
    with open(newfile, 'w') as f:
        f.write(new_descr.rstrip())

    oldfile = '/tmp/old_description_{}.txt'.format(uuid.uuid4())
    with open(oldfile, 'w') as f:
        f.write(old_descr)

    description_changed = bool(int(os.popen("diff "+oldfile+" "+newfile+" | wc -l").read()) > 0)
    os.remove(newfile)
    os.remove(oldfile)
    return description_changed



if __name__ == "__main__":
    main()
