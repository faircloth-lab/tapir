if (_WRAP_PASSTHROUGH_FLAG_) {
    fscanf (stdin, "String", treeString);
    if (_DO_TREE_REBALANCE_) {
        treeString = RerootTree (treeString, 0);
    }
    if (_KEEP_I_LABELS_) {
        INTERNAL_NODE_PREFIX = "intNode";
    }

    Tree givenTree = treeString;
    if (_KEEP_I_LABELS_) {
        INTERNAL_NODE_PREFIX = "Node";
    }
} else {

    if (!IS_TREE_PRESENT_IN_DATA) {

        fscanf (FILE_TREE, REWIND, "Raw", treeString);

        treeStringPattern = treeString$"^#NEXUS";
        if (treeStringPattern[0] >= 0) {
            ExecuteCommands (treeString);
            if (IS_TREE_PRESENT_IN_DATA == 0) {
                fprintf (stdout, "\nThis NEXUS file doesn't contain a valid tree block");
                return 1;
            }
            // always take the first tree found
            treeString = NEXUS_FILE_TREE_MATRIX[0][1];
        } else {
            treeStringPattern = treeString$"\\(";
            if (treeStringPattern[0] < 0) {
                fprintf (stdout, "\nThis doesn't seem to be a valid Newick string file. Can't find the opening parenthesis.\nHad:", treeString).
                return 1;
            } else {
                parenCounter = 1;
                strlength = Abs (treeString);
                cp = treeStringPattern[0] + 1;
                while ( cp < strlength && parenCounter ) {
                    cpc = treeString[cp];
                    if (cpc == "(") {
                        parenCounter = parenCounter + 1;
                    } else {
                        if (cpc == ")") {
                            parenCounter = parenCounter - 1;
                        }
                    }
                    cp = cp + 1;
                }

                if (parenCounter) {
                    fprintf (stdout, "\nThis doesn't seem to be a valid Newick string file. Can't match the parentheses.\nHad:", treeString).
                    return 1;
                }

                treeStringPattern = treeStringPattern[0];
                treeString = treeString[treeStringPattern][cp - 1];
            }
        }

        if (_DO_TREE_REBALANCE_) {
            treeString = RerootTree (treeString, 0);
        }

        if (_KEEP_I_LABELS_) {
            INTERNAL_NODE_PREFIX = "intNode";
        }

        Tree givenTree = treeString;
        if (_KEEP_I_LABELS_) {
            INTERNAL_NODE_PREFIX = "Node";
        }
    }
}
