#!/usr/bin/env python3
import sqlite3
import json
import flask
import re
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import re
import time
from plotly.validators.scatter.marker import SymbolValidator


import dash_table

from dash.exceptions import PreventUpdate

import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots

import pandas as pd
import numpy as np


##################################################################################################################
# loading files
##################################################################################################################

print("\n========= loading files =========\n")

## For each gene and year and cell type: how many papers?
gene_celltype_year_papers = pd.read_csv('data/gene_celltype_year_papers.csv')

## For each gene and year: how many papers?
gene_year_papers = pd.read_csv('data/gene_year_papers.csv')

## For each cell type & gene: What is the expression level and normalized number of papers?
gene_celltype_papercount_exp = pd.read_csv('data/gene_celltype_papercount_exp.csv')

## Genes shown initially when the page is loaded
starting_genes = open('data/starting_genes.csv').read().split(',')

## List of unique genes
all_genes = gene_celltype_year_papers['gene'].unique()


# create the gene dropdown for genes with known names (ie. name != ensembl_id)
# data_genelist is: {"ensg":"gene symbol"}
data_genelist = pd.read_csv('data/genelist.csv')
genes_dropdown_options = { v[0]: v[1] for v in data_genelist[ ['Ensembl.Gene.ID', 'Associated.Gene.Name'] ].values }
genes_dropdown_options_inv = {v: k for k, v in genes_dropdown_options.items()}

# add genes dictionary for genes display (this is a duplicate)   .... TODO, we should not have this twice ideally.
# Genes dict is: {"ensg":"gene symbol - long name"}
genes_dict = pd.read_csv('data/mus_musculus_genes_dictionary.csv')
genes_dict = { gene[1]: gene[2] for gene in genes_dict.values }
# genes_dict = { gene[1]: gene[2] + ' - ' + gene[3] for gene in genes_dict.values }

# for genename in data_genelist["Ensembl.Gene.ID"].tolist():
#     if genename not in genes_dict:
#         genes_dict[genename] = genename


for g in genes_dict:
    words = []
    lines = []
    for w in genes_dict[g].split(' '):
        words.append(w)
        if sum([len(word) for word in words]) > 40:
            lines.append(' '.join(words))
            words = []
    if len(words):
        lines.append(' '.join(words))
    genes_dict[g] = '<br>'.join(lines)




##################################################################################################################
# Helper functions
##################################################################################################################
# def parse_genes(genes_textbox):
#
#     '''
#     sanitizes the list of genes provided in the textbox
#     converting it into a python array of valid ensembl genes
#     '''
#
#     # split at commas and remove trailing spaces
#     genes = [ g.strip() for g in genes_textbox.split(',') ]
#
#     # replace gene names with their ensembl ids
#     genes = [ genes_dropdown_options[g] if g in genes_dropdown_options else g for g in genes ]
#
#     # removes mistyped and duplicate genes
#     genes = list(set([ g for g in genes if g in genes_dictionary ]))
#
#     return genes
###################################
def convert_genenames_to_ensg(genes):
    # replace gene names with their ensembl ids
    #genes = [ genes_dropdown_options[g] if g in genes_dropdown_options else g for g in genes ]
    genes = [ genes_dropdown_options_inv[g] if g in genes_dropdown_options_inv else g for g in genes ]
    # removes mistyped and duplicate genes
    #genes = list(set([ g for g in genes if g in genes_dictionary ]))
    genes = list(set([ g for g in genes if g in genes_dict ]))
    return genes

###################################
def convert_ensg_to_genenames(genes):
    genes = [ genes_dropdown_options[g] for g in genes ]
    return genes

###################################
###################################
def parse_genes(genes_textbox):
    '''
    sanitizes the list of genes provided in the textbox
    converting it into a python array of valid ensembl genes
    '''

    # split at commas and remove trailing spaces
    genes = [ g.strip() for g in genes_textbox.split(',') ]
    print("type"+type(genes))
    return convert_genenames_to_ensg(genes)

##################################################################################################################
# plotting functions
##################################################################################################################

def scatterplot(celltype = '', color = '', selected_genes = [], coord_data_plot = pd.DataFrame(), celltype_dependence = False):
    # print(type(selected_genes))
        
    t1 = time.time()

    if not len(selected_genes):
        fig = {}
        return fig
    else:

        selected_genes = list(selected_genes.split(","))
        selected_genes = convert_ensg_to_genenames(selected_genes)

        if celltype_dependence:
            xcol = "x_" + celltype
            ycol = "y_" + celltype
        else:
            xcol = "x"
            ycol = "y"


        print(selected_genes)
        coord_data_plot = coord_data_plot.dropna()

        raw_symbols = SymbolValidator().values
        marker_symbols = [raw_symbols[31] if v in selected_genes else raw_symbols[1] for i,v in enumerate(coord_data_plot["gene"].tolist())]

        marker_sizes = [45 if v in selected_genes else 5 for i,v in enumerate(coord_data_plot["gene"].tolist())]

        xaxis = coord_data_plot[xcol].values.tolist()
        yaxis = coord_data_plot[ycol].values.tolist()
        markercolor = coord_data_plot[color].values.tolist()
        textvalues = coord_data_plot["gene"].values.tolist()

        fig = go.Figure(go.Scatter(x = xaxis, y = yaxis, mode = "markers", 
            marker_color = markercolor, text = textvalues, opacity = 1.0,
        marker_symbol = marker_symbols, marker=dict(size=marker_sizes)))

        fig.update_layout( autosize= True)
        #fig.update_layout( autosize= False, width = 1800, height = 600)

        t2 = time.time()
        print(str(t2 -t1))
        return fig

##################################################################################################################
# starting server

print("\n========= starting server =========\n")

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']  #was not there before when nico made it. needed?
server = flask.Flask(__name__)
app = dash.Dash(
  __name__,
  server=server,
  routes_pathname_prefix='/bias/')#,
  # external_stylesheets=external_stylesheets)


app.config.suppress_callback_exceptions = True ###Note, dangerous -- see if needed. used if components are added after initialization
app.title = "Publicational bias visualization"

coordinate_list_data = pd.read_csv("data/list_coordinates.csv")
coordinate_list = {v[0]: v[1] for v in coordinate_list_data[["coord_id","coord_name"]].values}


conn = sqlite3.connect("data/totfeature.sqlite")
celltype_data = pd.read_sql_query("SELECT * from feature_matrix", conn)
conn.close()
celltype_list = celltype_data["ct"].unique()

feature_data = pd.read_csv("data/feature_long_name.csv")
feature_list = { v[0]: v[1] for v in feature_data[ ['feature_id', 'feature_long_name']].values }

# suggestions_names = genes_dropdown_options
# suggestions_names = convert_ensg_to_genenames(genes_dropdown_options)
# suggestions_ids = convert_genenames_to_ensg(genes_dropdown_options)
# suggestions_names.extend(suggestions_ids)

app.layout = html.Div([
    html.Div([
        html.Div([
            html.Div([
            html.Label('Gene selection'), # textbox for selecting genes using ensembl id or names; all plots updates are connected to this element
#            html.Datalist(
#                id='list-suggested-inputs',
#                children=[html.Option(value=word) for word in suggestions_names]),
             dcc.Input(
                id='gene-textbox',
                type='text',
                value='',
#               list='list-suggested-inputs',
                placeholder='Comma-separated list of genes to inspect',
                style={'width': '100%', 'height': '40px'}
            ),            # gene selection through dropdown; this will add the gene id to the textbox above
            html.Div([
                    dcc.Dropdown(
                    id='genes-dropdown',
                    value ='',
                    options=[ {'label': genes_dropdown_options[key], 'value':key} for key in genes_dropdown_options ],
                    placeholder='Select a gene using its name',)
            ], id='genes-dropdown-timestamp', n_clicks_timestamp = 0),],
            style={'width': '25%', 'display': 'inline-block'}),
            html.Div([
          #  html.Label('cell type'),
            html.Div([html.Label(['cell type'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='cell-type-selected',
                value= 'epithelial cell',
                options=[{'label': i, 'value': i} for i in celltype_list],
                placeholder = 'Cell type'),],
            style={'width': '25%', 'display': 'inline-block'}),
            html.Div([
            #html.Label('coordinate'),
            html.Div([html.Label(['coordinate'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='coord-id-selected',
                placeholder = 'coordinate',
                options=[{'label': i, 'value': i} for i in coordinate_list],
                value='coexp'),],
            style={'width': '25%', 'display': 'inline-block'}),
            html.Div([
            html.Div([html.Label(['color by'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='color-by-selected',
                placeholder = 'color by',
                options=[{'label': i, 'value': i} for i in feature_list],
                value='rank_pmid'),],
            style={'width': '25%', 'display': 'inline-block'}),
    ], style={'display': 'block', 'width': '100%'}),


    # without the following division the floating elements above will get overlapped by the plots
    # html.Div([], style={ 'padding': '10px 5px', 'display': 'inline-block', 'width': '100%' }),
        # hidden division with gene links to main databases;
        # upon gene identification will become visible and present the links to the user
        html.Div([
            html.Div([
                html.A(['Ensembl'], id='ensembl-link', href='', target='_blank')," | ",
                html.A(['UniProt'], id='uniprot-link', href='', target='_blank')," | ",
                html.A(['PubMed'], id='pubmed-link', href='', target='_blank')],),
                html.P([''])
        ], style={ 'width': '25%', 'display': 'block', 'float':'left'}),


        # umea university logo
        html.Div(
            [ html.Img(src='https://frontiersinblog.files.wordpress.com/2018/06/logo.png', style={
                     'height': '75%',
                     'float':'right',
                     'padding': '10px 50px'
                })
            ],
            style={
                'width': '25%',
                'height': '40px',
                'display': 'inline-block',
                'position': 'relative',
                'float':'right',
                'bottom': '0'
                }
        )

    ], style={
        'borderBottom': 'thin lightgrey solid',
        'backgroundColor': 'rgb(250, 250, 250)',
        'padding': '10px 5px',
        'width':'100%',
        'float':'left'
    }),
#
# without the following division the floating elements above will get overlapped by the plots
#  html.Div([], style={ 'padding': '10px 5px', 'display': 'inline-block', 'width': '100%' }),
        html.Div([
            dcc.Graph( id='scatter-plot', figure=scatterplot())
            ],#),
           style={'display': 'inline-block', 'width': '100%', 'height': '40%', 'margin': '0 auto',
          # style={'display': 'inline-block', 'width': '100%', 'height': '600px', 'margin': '0 auto',
           'padding': '10px 5px'}),
        ])



##################################################################################################################
# page callbacks
# test genes = ENSMUSG00000000126,ENSMUSG00000095309,ENSMUSG00000097090

##################################################################################################################
@app.callback(
    Output('gene-links-div', 'style'),
    [Input('gene-textbox', 'value')])
def update_gene_links_div(selected_genes):

    '''
    hide genes links if more than a gene is provided
    '''

    selected_genes = parse_genes(selected_genes)
    # print("Here"+type(selected_genes))

    if len(selected_genes) == 1:
        gene_links = {'display': 'block', 'padding':'30px 0 0 30px'}
    else:
        gene_links = {'display': 'none', 'padding':'30px 0 0 30px'}
    #
    # for i,v in enumerate(genes_dict):
    #     if i > 100:
    #         break
    #     else:
    #         print(i,v, genes_dict[v])
    return gene_links

#############################

################################################################################
@app.callback(
    Output('gene-textbox', 'value'),
    [Input('genes-dropdown', 'value')])
def update_genes_dropdown(dropdown_value):

    '''
    the genes textbox can be updated manually or using the dropdown
    '''

    if dropdown_value is None:
        return ''
    else:
        return dropdown_value
###################################################
@app.callback(Output('scatter-plot', 'figure'),
    [Input('gene-textbox', 'value'),
     Input('cell-type-selected', 'value'),
     Input('coord-id-selected', 'value'),
     Input('color-by-selected', 'value')])

def update_graph(selected_genes, celltype,coordid,color):

    t1 = time.time()
    conn = sqlite3.connect("data/coord_" + coordid + ".sqlite")
    coord_data = pd.read_sql_query("SELECT * from coord", conn)
    conn.close()
    
    conn = sqlite3.connect("data/totfeature.sqlite")
    celltype_data_view = pd.read_sql_query("SELECT * from feature_matrix where ct == \"" + celltype + "\"", conn)
    conn.close()

    #print(celltype_data_view)

    # print(coord_data.head(5))
    # selected_genes = parse_genes(selected_genes)
    # rows2view = [i for i,v in enumerate(celltype_data["ct"].tolist()) if v == celltype]
    # celltype_data_view = celltype_data.iloc[rows2view,:].copy()
    coord_data = coord_data.merge(celltype_data_view,left_on = "gene", right_on = "gene", indicator = True)

    #print(coord_data)

    celltype_dependence_data = pd.read_csv("data/list_coordinates.csv")
    celltype_dependence_data.index = celltype_dependence_data["coord_id"].tolist()
    celltype_dependence = celltype_dependence_data.loc[coordid,"perct"]


    if celltype_dependence:
        coord_data_plot = coord_data.loc[:,["gene","x_"+ celltype, "y_"+ celltype, color]].copy()
    else:
        coord_data_plot = coord_data

    print(coord_data_plot)

    return scatterplot(celltype, color, selected_genes, coord_data_plot, celltype_dependence)
#################################################################################
@app.callback(
    [Output('ensembl-link', 'href'),
    Output('uniprot-link', 'href'),
    Output('pubmed-link', 'href')],
    [Input('gene-textbox', 'value')])
def update_gene_links(gene):

    '''
    update gene links
    '''

    gene = gene.strip()

    links = (
        'https://www.ensembl.org/Mus_musculus/Gene/Summary?g={}'.format(gene),
        'https://www.uniprot.org/uniprot/?query={}&sort=score'.format(gene),
        'https://www.ncbi.nlm.nih.gov/search/all/?term={}'.format(gene),
    )

    return links
###################################



# run the server

# run the app on "python app.py";
# default port: 8050
if __name__ == '__main__':
    app.run_server(debug = True)

app = dash.Dash(__name__)
#viewer.show(app)
