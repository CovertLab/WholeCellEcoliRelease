#! /usr/bin/env bash
# Remove compiled python files, run make compile and check requirements
# to maintain proper environment for wcm.

echo -e "\nRemoving *.pyc..."
find . -not \( -path ./out -prune \) -not \( -path ./.git -prune \) -name "*.pyc" -exec rm {} \;

echo "Running make compile..."
make compile

if [ -e requirements.txt ]; then
    echo -e "\nRequirements diff (requirements.txt vs current pips):"
    diff --ignore-case <(sed 's/ *#.*//;s/^ *--.*//;/^$/d' requirements.txt | sort --ignore-case) \
      <(pip freeze 2>/dev/null | sort --ignore-case) -yB --suppress-common-lines
fi

# Exit normally so git rebase works even if diff finds a diff (error exit code)
exit 0
