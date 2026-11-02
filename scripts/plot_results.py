"""Export an English comparison figure; all distances are normalized by GH."""
from pathlib import Path
import argparse,csv
ROOT=Path(__file__).resolve().parents[1]

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results/figures')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    with (ROOT/'results/unfinished_5_comparison.csv').open() as stream:rows=list(csv.DictReader(stream))
    print('TABLE_COLUMNS',','.join(rows[0]),flush=True)
    fields=dict(lower='lambda2_certified_lower_log2',candidate='projected_value_log2',paper='revised_pdf_lambda2_upper_log2',bkz='local_BKZ_lambda2_upper_log2')
    # Column names are checked against the table schema before plotting.
    if any(v not in rows[0] for v in fields.values()):raise ValueError('Plot/table schema mismatch')
    fig,ax=plt.subplots(figsize=(8,3.8),layout='constrained')
    for i,row in enumerate(rows):
        gh=float(row['prediction_GH_log2']);lo=2**(float(row[fields['lower']])-gh);hi=2**(float(row[fields['paper']])-gh)
        ax.plot([lo,hi],[i,i],color='0.6',lw=3)
        ax.scatter(lo,i,marker='|',s=150,color='tab:blue',label='Certified lower bound' if i==0 else None)
        ax.scatter(hi,i,marker='|',s=150,color='tab:red',label='Revised paper upper (rounded)' if i==0 else None)
        ax.scatter(2**(float(row[fields['candidate']])-gh),i,marker='o',s=32,color='tab:green',label='Projected candidate (not exact)' if i==0 else None)
        if row[fields['bkz']]:ax.scatter(2**(float(row[fields['bkz']])-gh),i,marker='x',s=35,color='tab:purple',label='Original-lattice BKZ upper' if i==0 else None)
    ax.axvline(1,color='black',ls='--',lw=1,label='GH prediction')
    ax.set_yticks(range(len(rows)),[r['case'] for r in rows]);ax.invert_yaxis()
    ax.set_xlabel('Length / GH prediction');ax.legend(fontsize=7,loc='upper left',bbox_to_anchor=(1,1))
    fig.savefig(args.output/'unfinished_bounds.png',dpi=180,metadata={'Software':'Matplotlib'})
    fig.savefig(args.output/'unfinished_bounds.pdf',metadata={'Author':'','Creator':'Matplotlib','CreationDate':None,'ModDate':None})
    plt.close(fig)

if __name__=='__main__':main()
