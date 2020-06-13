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
import math

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

## Load coordinate system names
coordinate_list_data = pd.read_csv("data/list_coordinates.csv")
coordinate_list = {v[0]: v[1] for v in coordinate_list_data[["coord_id","coord_name"]].values}

## Load cell types
conn = sqlite3.connect("file:data/totfeature.sqlite?mode=ro", uri=True)
celltype_data = pd.read_sql_query("SELECT DISTINCT ct from feature_matrix ORDER BY ct", conn) 
celltype_list = celltype_data["ct"]
conn.close()

## Load features (to be shown as color)
feature_data = pd.read_csv("data/feature_long_name.csv")
feature_list = { v[0]: v[1] for v in feature_data[ ['feature_id', 'feature_long_name']].values }

## Load cell types
celltype_dependence_data = pd.read_csv("data/list_coordinates.csv")
celltype_dependence_data.index = celltype_dependence_data["coord_id"].tolist()


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
    #print("type"+type(genes))
    return convert_genenames_to_ensg(genes)







##################################################################################################################
# Function: Make the histogram of #citations over time for a given gene
##################################################################################################################
def histogram_citationsperyear(gene_id):
  
    ## Load histogram of citations for this gene
    conn = sqlite3.connect("file:data/citations_per_year.sqlite?mode=ro", uri=True)
    citationsperyear = pd.read_sql_query("SELECT * from citationsperyear where ensembl == ?", conn, params=(gene_id,))
    conn.close()
        
    ##Check if there is any data to plot
    if citationsperyear.shape[0]==0:
        fig = go.Figure()
        fig.update_layout( autosize= False, width = 400, height = 100, margin={'t':0, 'b':0,'l':0, 'r':0})
        return fig

    ##Check if any genes should be highlighted
    fig = go.Figure()
    fig.add_trace(go.Scatter(
      x=citationsperyear["year"].tolist(), 
      y=citationsperyear["citations"].tolist(), 
      fill='tozeroy')) # fill down to xaxis
        
    fig.update_layout( autosize= False, width = 400, height = 100, margin={'t':0, 'b':0,'l':0, 'r':0})
    
    return fig



##################################################################################################################
##### The main window layout
##################################################################################################################

print("\n========= starting server =========\n")

server = flask.Flask(__name__)
app = dash.Dash(
  __name__,
  server=server,
  routes_pathname_prefix='/genepub/')

app.config.suppress_callback_exceptions = True ###Note, dangerous -- see if needed. used if components are added after initialization
app.title = "Data viewer: 10 reasons to study a gene"
app.layout = html.Div([
    html.Div([
        html.Div([
    
            # textbox for selecting genes using ensembl id or names; all plots updates are connected to this element
            html.Label('Selected gene:'), 
            dcc.Input(
                id='gene-textbox',
                type='text',
                value='',
                #list='list-suggested-inputs',  #Don't suggest EnsemblIDs
                placeholder='Comma-separated list of genes to inspect',
                style={'width': '100%', 'height': '40px'}
            ),            
            # gene selection through dropdown; this will add the gene id to the textbox above
            html.Div([
                    dcc.Dropdown(
                    id='genes-dropdown',
                    value ='',
                    options=[ {'label': genes_dropdown_options[key], 'value':key} for key in genes_dropdown_options ],
                    placeholder='Select a gene using its name',)
            ], id='genes-dropdown-timestamp', n_clicks_timestamp = 0),
            
            
            html.Div([html.Label(['Cell type:'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='cell-type-selected',
                value= 'epithelial cell',
                options=[{'label': i, 'value': i} for i in celltype_list],
                placeholder = 'Cell type'),
                
                
            html.Div([html.Label(['Coordinates:'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='coord-id-selected',
                placeholder = 'coordinate',
                options=[{'label': i, 'value': i} for i in coordinate_list],
                value='coexp'),
                
                
            html.Div([html.Label(['Color by:'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='color-by-selected',
                placeholder = 'color by',
                options=[{'label': i, 'value': i} for i in feature_list],
                value='rank_pmid'),
    
            html.Div([
                html.A([
                    html.Img(src=app.get_asset_url('MIMS_logo_blue.svg'), style={
                           'height': '30px',
                           #'float':'right',
                           'padding': '10px 10px'}),
                ], href='http://www.mims.umu.se/')
            ], style={'text-align':'right'})
                  
        ], style={
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 5px',
            'width':'100%'
        }),
        
        #html.Br(),

        #####################################################
        ########## Gene information panel ###################
        #####################################################
        html.Div([
          
            #html.H1(html.Label(id='geneinfo-symbol')),
            html.H4(html.Label(id='geneinfo-longname')),
            html.Label(id="geneinfo-firstcited"),
            html.Label(id="geneinfo-numcitations"),

            html.Div([
              dcc.Graph( id='citationsperyear-histogram')
            ],style={
              'display': 'inline-block',
              'margin': 'auto'
            }),
            
            html.Div([
                html.A(['Ensembl'],    id='ensembl-link',   href='', target='_blank')," | ",
                html.A(['UniProt'],    id='uniprot-link',   href='', target='_blank')," | ",
                html.A(['PubMed'],     id='pubmed-link',    href='', target='_blank')," | ",
                html.A(['Genecards'],  id='genecards-link', href='', target='_blank')
            ])


        ], style={
            'margin-top':'50px',
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 5px',
            'width':'100%'
        }, id='geneinfo-div'),
    ], style={'float':'right','width':'25%', 'padding':'20px'}),
    html.Div([
        dcc.Graph( id='scatter-plot')#, figure=scatterplot())
    ],style={
      'display': 'inline-block',
      'margin': '0 auto'
    }),
],style={
  'position': 'inline-block', 
  'width': '95%', 
  'height': '95%', 
  'margin': '0 auto', #supposed to center it... but meah.
  'padding':'0',
  'overflow':'hidden'
}) 


##################################################################################################################
# page callbacks
# test genes = ENSMUSG00000000126,ENSMUSG00000095309,ENSMUSG00000097090
##################################################################################################################




################################################################################
@app.callback(
    Output('gene-textbox', 'value'),
    [Input('genes-dropdown', 'value')])
def update_genes_dropdown(dropdown_value):
    ##the genes textbox can be updated manually or using the dropdown
    if dropdown_value is None:
        return ''
    else:
        return dropdown_value
        
        
##################################################################################################################
# Function: Make the scatter plot for all the genes
##################################################################################################################
@app.callback(Output('scatter-plot', 'figure'),
    [Input('gene-textbox', 'value'),
     Input('cell-type-selected', 'value'),
     Input('coord-id-selected', 'value'),
     Input('color-by-selected', 'value')])
def update_graph(selected_genes, celltype,coordid,color):

    ##Check if the coordinates are per cell type; otherwise just use the columns x,y
    celltype_dependence = celltype_dependence_data.loc[coordid,"perct"]
    if celltype_dependence:
        xcoord_name = "x_"+celltype
        ycoord_name = "y_"+celltype
    else:
        xcoord_name = "x"
        ycoord_name = "y"

    ##Load coordinates
    conn = sqlite3.connect("file:data/coord_" + coordid + ".sqlite?mode=ro", uri=True)
    coord_data = pd.read_sql_query("SELECT gene, `"+xcoord_name+"` as x, `"+ycoord_name+"` as y from coord", conn)
    conn.close()

    ##Load the feature we need
    conn = sqlite3.connect("file:data/totfeature.sqlite?mode=ro", uri=True)
    feature_data = pd.read_sql_query("SELECT gene, `"+color+"` as feature from feature_matrix where ct == ?", conn, params=(celltype,))
    conn.close()

    ##Merge coordinates and features
    coord_data_plot = coord_data.merge(feature_data,left_on = "gene", right_on = "gene", indicator = True)

    ##Check if there is any data to plot, otherwise return an empty graph
    if coord_data_plot.shape[0]==0:
        fig = go.Figure()
        fig.update_layout( autosize= False, width = 800, height = 800)
        return fig

    ##Check if any genes are selected and should be highlighted
    if not len(selected_genes):
        selected_genes = []
    else:
        selected_genes = list(selected_genes.split(","))
        selected_genes = convert_ensg_to_genenames(selected_genes)

    ##Drop crappy values - usually missing features
    coord_data_plot = coord_data_plot.dropna()

    #Extract X,Y coordinates of selected genes
    vlines = [coord_data_plot["x"].values[i] for i,v in enumerate(coord_data_plot["gene"].tolist()) if v in selected_genes]
    hlines = [coord_data_plot["y"].values[i] for i,v in enumerate(coord_data_plot["gene"].tolist()) if v in selected_genes]

    #Extarct X,Y coordinates of all genes, their feature value (color), and the name
    xaxis = coord_data_plot["x"].values.tolist()
    yaxis = coord_data_plot["y"].values.tolist()
    markercolor = coord_data_plot["feature"].values.tolist()
    textvalues = coord_data_plot["gene"].values.tolist()

    #Create the basic plot
    fig = go.Figure(
        go.Scatter(
            x = xaxis, 
            y = yaxis, 
            mode = "markers",
            marker_color = markercolor, 
            text = textvalues,
            opacity = 1.0))

    #Add cross-hairs to all the selected genes
    shapes_y=[{'type': 'line',
                        'y0':y_intercept,
                        'y1':y_intercept,
                        'x0':str(min(xaxis)), 
                        'x1':str(max(xaxis)),
                        'line': {'color': 'black', 'width': 1, 'dash': 'dot'}} 
                   for i, y_intercept in enumerate(hlines)]
    shapes_x=[{'type': 'line',
                        'x0':x_intercept,
                        'x1':x_intercept,
                        'y0':str(min(yaxis)), 
                        'y1':str(max(yaxis)),
                        'line': {'color': 'black', 'width': 1, 'dash': 'dot'}} 
                   for i, x_intercept in enumerate(vlines)]
    fig.layout.update(shapes=shapes_x+shapes_y)
        
    fig.update_layout( autosize= False, width = 800, height = 800, margin={'t':0, 'b':0,'l':0, 'r':0})
    return fig    



##################################################################################################################
##### Callback: Update gene information box ... links
##################################################################################################################
@app.callback(
    [
    Output('geneinfo-div',          'style'),
    #Output('geneinfo-symbol',       'children'),
    Output('geneinfo-longname',     'children'),
    Output('geneinfo-firstcited',   'children'),
    Output('geneinfo-numcitations', 'children'),
    
    Output('citationsperyear-histogram', 'figure'),

    Output('ensembl-link',   'href'),
    Output('uniprot-link',   'href'),
    Output('pubmed-link',    'href'),
    Output('genecards-link', 'href')
    ],
    [Input('gene-textbox', 'value')])  
def update_gene_links(gene):

    ##hide genes links if more than a gene is provided
    selected_genes = parse_genes(gene)
    gene_id = selected_genes
    gene_symbol = convert_ensg_to_genenames(selected_genes)
    geneinfo_longname=""
    geneinfo_firstcited=""
    geneinfo_numcitations=""
    retfig_citationsperyear=""

    if len(selected_genes) == 1:
        gene_symbol = gene_symbol[0]
        gene_id = gene_id[0]

        ### Make geneinfo box visible
        style= {
            'margin-top':'50px',
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 5px',
            'width':'100%'
        }

        ### Pull out information about the gene
        conn = sqlite3.connect("file:data/geneinfo.sqlite?mode=ro", uri=True)
        geneinfo_data = pd.read_sql_query("SELECT * from geneinfo WHERE ensembl=?", conn, params=(gene_id,))
        conn.close()
        if geneinfo_data.shape[0]>0:
            geneinfo_data = geneinfo_data.to_dict()
            geneinfo_longname     = geneinfo_data["description"][0]
            
            if geneinfo_data["firstyear"][0]<0:
                geneinfo_firstcited   = "First cited: N/A"
            else:
                geneinfo_firstcited   = "First cited: "+str(int(geneinfo_data["firstyear"][0]))
            
            geneinfo_numcitations = "#Citations: "+str(geneinfo_data["numcitations"][0])

        retfig_citationsperyear = histogram_citationsperyear(gene_id)

        #print(citationsperyear)


    else:
        gene_symbol=""
        style= {'display': 'none', 'padding':'30px 0 0 30px'}


    ##update gene links
    gene = gene.strip()
    info = (
        style,
        
        #gene_symbol,
        geneinfo_longname,
        geneinfo_firstcited,
        geneinfo_numcitations,
        
        retfig_citationsperyear,

        'https://www.ensembl.org/Mus_musculus/Gene/Summary?g={}'.format(gene),
        'https://www.uniprot.org/uniprot/?query={}&sort=score'.format(gene),
        'https://www.ncbi.nlm.nih.gov/search/all/?term={}'.format(gene),
        'https://www.genecards.org/cgi-bin/carddisp.pl?gene={}&keywords={}'.format(gene_symbol,gene_symbol)
    )

    return info


###################################

# run the app on "python app.py";
# default port: 8050
if __name__ == '__main__':
    app.run_server(debug = True)

app = dash.Dash(__name__)
#viewer.show(app)














