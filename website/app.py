#!/usr/bin/env python3
import sqlite3
import json
import flask
import re
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import time

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
genes_dict = { gene[1]: gene[2] + ' - ' + gene[3] for gene in genes_dict.values }
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
    return convert_genenames_to_ensg(genes)

##################################################################################################################
# plotting functions
##################################################################################################################


#######################################################
## Heatmap plot: Gene paper count, year vs cell type ##
#######################################################
def plot_gene_celltype_papers(gene):

    '''
    heatmap plot of paper counts per gene and celltype;
    x: celltype, y: year, year will be plotted from min_year (exclusive) to max_year (inclusive);
    the provided gene should have already been sanitized
    '''

    if len(gene)==0:
        return {}
        # return {'display': 'none'}
    #{'display': 'block', 'width': '100%'}

    gene = convert_ensg_to_genenames(gene)[0]

    ## Get the data from the sqlite database
    conn = sqlite3.connect("data/gene_celltype_papercount_year.sqlite")
    df = pd.read_sql_query("SELECT * from "+gene, conn)
    conn.close()
    df = df.pivot(index='Var2', columns='Var1', values='value')
    years = df.columns.values
    celltypes = df.index.values
    expression_matrix = df.to_numpy()

    # generate heatmap using the computed expression matrix
    fig = go.Figure(
        layout = {
            'title': gene,
            'height': len(celltypes)*20 + 200, # plot height should be proportional to the number of entries
            'yaxis': {
                'side': 'right',
                'autorange': 'reversed', # put first celltypes on top rather than on the bottom
            },
        },
        data = go.Heatmap(
            z = expression_matrix,
            x = years,
            y = celltypes,
            colorbar = {'x': -.1},
            showscale = False,
        )
    )

    # fig.update_yaxes()

    return fig



###########################################
## Plot: Gene expression vs paper count  ##
###########################################
def plot_paper_gene_exp(genelist):
    print(genelist, convert_ensg_to_genenames(genelist))
    if len(genelist)==0:
        return {}
        # return {'display': 'none'}

    conn = sqlite3.connect("data/gene_celltype_papercount.sqlite")
    df = pd.read_sql_query("SELECT * from gene_cell_count", conn)
    conn.close()

    genes = [convert_ensg_to_genenames(genelist)[0]]
    subplot_titles = [genename for igene, genename in enumerate(genes)]
    fig = make_subplots(rows=len(genelist), cols=1, subplot_titles = tuple(subplot_titles))
    for i, gene in enumerate(genes):
        gene_mask = [ig for ig, genename in enumerate(df['gene'].values) if genename == gene]
        dfplot = df.iloc[gene_mask,:].copy()

        xaxis = dfplot["term"].values.tolist()
        yaxis = dfplot["papercount"].values.tolist()

        fig.add_trace(go.Scatter(x=xaxis, y=yaxis))



    return fig


##################################################################################################################
# starting server

print("\n========= starting server =========\n")

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']  #was not there before when nico made it. needed?
server = flask.Flask(__name__)
app = dash.Dash(
  __name__,
  server=server,
  routes_pathname_prefix='/bias/',
  external_stylesheets=external_stylesheets)


app.config.suppress_callback_exceptions = True ###Note, dangerous -- see if needed. used if components are added after initialization
app.title = "Publicational bias visualization"


##################################################################################################################
# page layout

app.layout = html.Div([
    html.Div([
        html.Div([

            html.Label('Gene selection'),

            # textbox for selecting genes using ensembl id or names;
            # all plots updates are connected to this element
            dcc.Input(
                id='gene-textbox',
                type='text',
                value='',
                placeholder='Comma-separated list of genes to inspect',
                style={'width': '100%'}
            ),

            # gene selection through dropdown; this will add the gene id to the textbox above
            html.Div([
                dcc.Dropdown(
                    id='genes-dropdown',
                    options=[ {'label': genes_dropdown_options[key], 'value':key} for key in genes_dropdown_options ],
                    placeholder='Select a gene using its name',
                )
            ], id='genes-dropdown-timestamp', n_clicks_timestamp = 0),

        ],
        style={'width': '75%', 'display': 'inline-block', 'float':'left'}),

        # umea university logo
        html.Div(
            [ html.Img(src='https://frontiersinblog.files.wordpress.com/2018/06/logo.png', style={
                     'height': '100%',
                     'float':'right',
                     'padding': '10px 30px'
                })
            ],
            style={
                'width': '25%',
                'height': '40px',
                'display': 'inline-block',
                'position': 'relative',
                'float':'left',
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

    # without the following division the floating elements above will get overlapped by the plots
    html.Div([], style={ 'padding': '10px 5px', 'display': 'inline-block', 'width': '100%' }),

    html.Div([
        html.Div([
            #Plot heatmap for paper count, celltype vs year
            dcc.Graph( id='plot-pc-year-celltype', figure=plot_gene_celltype_papers(convert_genenames_to_ensg(["Gata4"])) )
        ], style={'display': 'block'}),
    ], style={ 'padding': '10px 5px', 'width': '100%', 'display': 'inline-block' }),

    html.Div([

        # genes plots; initialize as hidden
        html.Div([
		#plot number of paper vs gene expression
            #dcc.Graph( id='genes-plot', figure=genes_scatterplots() )
            dcc.Graph( id='plot-paper-gene-exp', figure=plot_paper_gene_exp(convert_genenames_to_ensg(["Gata4"])))
        ], style={'width': '100%'}, id='genes-plot-div'),
    ], style={ 'padding': '10px 5px', 'width': '100%', 'display': 'inline-block'}),

],style={'max-width': '1200px', 'margin': '0 auto'})


##################################################################################################################
# page callbacks
##################################################################################################################






###################################
###################################
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


###################################
###################################
@app.callback(
    [Output('plot-pc-year-celltype', 'figure')],
#    [Output('genes-plot', 'figure'),
#    Output('genes-plot-div', 'style')],
    [Input('gene-textbox', 'value')])
def update_plots(selected_genes):

    '''
    plots update, triggered by genes selection
    '''

    selected_genes = parse_genes(selected_genes)
    print(selected_genes)

    plot_pc_year_celltype = plot_gene_celltype_papers(selected_genes)

    #return umap_plot_cells(selected_genes), genes_scatterplots(selected_genes), genes_plot
    #return genes_scatterplots(selected_genes), genes_plot
    return [plot_pc_year_celltype]

###################################
###################################
@app.callback(
    [Output('plot-paper-gene-exp', 'figure')],
    [Input('gene-textbox', 'value')])
def update_plots(selected_genes):

    '''
    plots update, triggered by genes selection
    '''

    selected_genes = parse_genes(selected_genes)
    #print(selected_genes)

    plotpaper_gene_exp = plot_paper_gene_exp(selected_genes)

    return [plotpaper_gene_exp]

##################################################################################################################
# run the server

# run the app on "python app.py";
# default port: 8050
if __name__ == '__main__':
    app.run_server(debug = True)
