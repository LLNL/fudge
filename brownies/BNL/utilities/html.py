import collections

from brownies.BNL.inter import datatables


# -------------------------------------------------------------------------------
# Simple HTML helpers
# -------------------------------------------------------------------------------
def attributeFixer(k):
    k = str(k)
    if k.startswith('_'):
        k = k[1:]
    return k


def XML(tag, content, **kwds):
    result = ' '.join(["<" + str(tag)] + [attributeFixer(k) + '="' + str(v) + '"' for k, v in kwds.items()])
    if content is None:
        return result + '/>'
    return result + ">" + str(content) + "</" + str(tag) + ">"


def H1(x, **kwds):
    return XML('H1', x, **kwds)


def H2(x, **kwds):
    return XML('H2', x, **kwds)


def H3(x, **kwds):
    return XML('H3', x, **kwds)


def H4(x, **kwds):
    return XML('H4', x, **kwds)


def SCRIPT(x, **kwds):
    return XML('SCRIPT', x, **kwds)


def OBJECT2HTML(d, nametag="<strong>%s: </strong>", itemtag='<li>%s</li>', valuetag="%s", blocktag=('<ul>', '</ul>')):
    def islistlike(x):
        return isinstance(x, list) or isinstance(x, tuple) or isinstance(x, set)

    def isdictlike(x):
        return isinstance(x, dict) or isinstance(x, collections.OrderedDict)

    def dorecursion(x):
        return islistlike(x) or isdictlike(x)

    r = []
    if isdictlike(d):
        r.append(blocktag[0])
        for k, v in d.items():
            name = nametag % k
            if dorecursion(v):
                r.append(itemtag % name)
                OBJECT2HTML(r, v)
            else:
                value = valuetag % v
                r.append(itemtag % (name + value))
        r.append(blocktag[1])
    elif islistlike(d):
        r.append(blocktag[0])
        for i in d:
            if dorecursion(i):
                r.append(itemtag % " - ")
                OBJECT2HTML(r, i)
            else:
                r.append(itemtag % i)
        r.append(blocktag[1])
    elif isinstance(d, datatables.DataTable):
        r.append(d.to_html())
    return r


__symbolDict = {x[8:]: x[0:7] for x in '''&#0913; Alpha
&#0914; Beta
&#0915; Gamma
&#0916; Delta
&#0917; Epsilon
&#0918; Zeta
&#0919; Eta
&#0920; Theta
&#0921; Iota
&#0922; Kappa
&#0923; Lambda
&#0924; Mu
&#0925; Nu
&#0926; Xi
&#0927; Omicron
&#0928; Pi
&#0929; Rho
&#0931; Sigma
&#0932; Tau
&#0933; Upsilon
&#0934; Phi
&#0935; Chi
&#0936; Psi
&#0937; Omega
&#0945; alpha
&#0946; beta
&#0947; gamma
&#0948; delta
&#0949; epsilon
&#0950; zeta
&#0951; eta
&#0952; theta
&#0953; iota
&#0954; kappa
&#0955; lambda
&#0956; mu
&#0957; nu
&#0958; xi
&#0959; omicron
&#0960; pi
&#0961; rho
&#0962; final sigma
&#0963; sigma
&#0964; tau
&#0965; upsilon
&#0966; phi
&#0967; chi
&#0968; psi
&#0969; omega
&#0976; beta symbol
&#0977; theta symbol
&#0978; upsilon with hook symbol
&#0981; phi symbol
&#0982; pi symbol
&#0986; Stigma
&#0987; stigma
&#0988; Digamma
&#0989; digamma
&#1008; kappa symbol
&#1009; rho symbol
&#1012; Theta symbol
&#1013; lunate epsilon symbol'''.split('\n')}


def SYMBOL(key):
    return __symbolDict[key]


def TABLE(tab, row_labels=None):
    def TH(h):
        if hasattr(h, 'unit') and h.unit == '':
            return '<TH>%s</TH>' % (str(h.name).replace('<', '&lt;').replace('>', '&gt;'))
        return '<TH>%s</TH>' % (str(h).replace('<', '&lt;').replace('>', '&gt;'))

    TRs = []
    if row_labels and len(tab.data) != len(row_labels):
        raise ValueError("Number of rows in table and row_labels not equal")
    if row_labels:
        TRs.append('<TR>%s</TR>' % ' '.join([TH(x) for x in [''] + tab.columns]))
    else:
        TRs.append('<TR>%s</TR>' % ' '.join([TH(x) for x in tab.columns]))
    for irow, row in enumerate(tab):
        if row_labels:
            TR = ['<TD>%s</TD>' % str(x) for x in [row_labels[irow]] + row]
        else:
            TR = ['<TD>%s</TD>' % str(x) for x in row]
        TRs.append('<TR>%s</TR>' % ' '.join(TR))
    return "<TABLE>%s</TABLE>" % '\n'.join(TRs)


def IMG(src):
    return '<IMG src="%s">' % src

