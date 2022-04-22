from pathlib import Path
import csv

def get_rows(db_path, params, size=2):
    res = list()
    with open(db_path, 'rt', newline='') as db_file:
        reader = csv.DictReader(db_file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        for row in reader:
            if int(row['size']) == size:
                res.append([row[param] for param in params])
    return res




























'''

    clique db query system

    users should be able to easily and concisely make complex queries to db,
    since several fields require string parsing / iteration beyond simply comparing values,
    building a query function which takes a query using a more expressive and job-appropriate
    query language and parses it to automatically check the query predicates given against every row in the db
    is the most elegant solution

    fields: id,size,clique,resid,oldresid,layerinfo,pdbname




    -size 1;2;3 / +size 1;2;-3 / clique (1;3;4); -(1;3)

    query language grammar:

    Expr ->
        Clause '/' Expr
        | Clause
    Clause ->
        Pos_Clause
        | Neg_Clause
    Pos_Clause ->
        Filter
        | '+' Filter
    Neg_Clause ->
        '-' Filter
    Filter ->
        Single_Filter
        | Multi_Filter
    Single_Filter ->
        'size' Single_Rules
        | 'id' Single_Rules
        | 'pdbname' Single_Rules
    Multi_Filter ->
        'clique' Multi_Rules
        | 'resid' Multi_Rules
        | 'oldresid' Multi_Rules
        | 'layerinfo' Multi_Rules
    Single_Rules ->
        Single_Rule ';' Single_Rules
        | Single_Rule
    Multi_Rules ->
        Multi_Rule ';' Multi_Rules
        | Multi_Rule
    Single_Rule ->
        Pos_Single_Rule
        | Neg_Single_Rule
    Multi_Rule ->
        Pos_Multi_Rule
        | Neg_Multi_Rule
    Pos_Single_Rule ->
        Value
        | '+' Value
    Neg_Single_Rule ->
        '-' Value
    Pos_Multi_Rule ->
        Multi_Value
        | '+' Multi_Value
    Neg_Multi_Rule ->
        '-' Multi_Value
    Multi_Value ->
        Value
        | Value_Set
    Value_Set ->
        '(' Inner_Value_Set ')'
    Inner_Value_Set ->
        Value ',' Inner_Value_Set
        | Value
    Value ->
        INT
        | STRING
        | IntRange
    IntRange ->
        INT ':' INT

    fields: id,size,clique,resid,oldresid,layerinfo,pdbname

    query language syntax grammar:

    Query ->
        Clause '/' Query
        | Clause
    Clause ->
        '+' Sub_Clause
        | '-' Sub_Clause
        | Sub_Clause
    Sub_Clause ->
        'id' Rules
        | 'size' Rules
        | 'clique' Rules
        | 'resid' Rules
        | 'oldresid' Rules
        | 'layerinfo' Rules
        | 'pdbname' PDBNAME_Rules
    PDBNAME_Rules ->
        REGEX_STRING ';' PDBNAME_Rules
        | REGEX_STRING
    Rules ->
        Rule ';' Rules
        | Rule
    Rule ->
        '+' Sub_Rule
        | '-' Sub_Rule
        | Sub_Rule
    Sub_Rule ->
        '(' Inner_Value_Set ')'
        | INT ':' INT
        | INT
        | RESIDUE
    Inner_Value_Set ->
        Inner_Value_Set_INT
        | Inner_Value_Set_RES
    Inner_Value_Set_INT ->
        INT ',' Inner_Value_Set_INT
        | INT
    Inner_Value_Set_RES ->
        RESIDUE ',' Inner_Value_Set_RES
        | RESIDUE

    ideas for query syntax tokenizer:
        - hardcoded tokenizer, almost done,
            responds poorly to changes in
            syntax grammar
        - general top-down matcher / parser,
            less code, adapts to changes in
            grammar, most of it can be reused for
            parsing tokens against language grammar

'''




