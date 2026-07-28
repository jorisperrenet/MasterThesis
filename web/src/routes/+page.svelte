<script lang="ts">
  import { base } from '$app/paths';
  import SiteMeta from '$lib/site-kit/SiteMeta.svelte';

  const title = 'Decoding CSIDH: Isogeny-Based Cryptography Explained — Joris Perrenet';
  const description = 'A self-contained guide to CSIDH, elliptic curves, isogenies, ideals and worked SageMath examples.';
  const canonical = 'https://jorisperrenet.com/MasterThesis/';
  const pdfUrl = `${canonical}decoding-csidh.pdf`;
  const abstract = `In today's digital age, private communication is anchored in robust encryption techniques. As the horizon of quantum computing draws nearer, current encryption methods may soon become vulnerable to quantum attacks, pointing to the need for post-quantum cryptography. The CSIDH encryption scheme, grounded in isogeny-based cryptography, emerges as a compelling candidate for such post-quantum cryptography. Understanding CSIDH requires a solid grasp of finite fields, number fields, and particularly elliptic curves. This thesis strives to break down these topics to make the innovative CSIDH approach more accessible, encouraging further exploration in the field of isogeny-based cryptography to ensure digital privacy for years to come.`;
  const bibtex = `@mastersthesis{perrenet2024decoding,
  author = {Perrenet, L. J.},
  title  = {Decoding CSIDH: A Guide to Isogeny-Based Cryptography},
  school = {Leiden University},
  year   = {2024},
  month  = jul,
  url    = {${pdfUrl}}
}`;
  let citationCopied = $state(false);

  async function copyCitation() {
    await navigator.clipboard.writeText(bibtex);
    citationCopied = true;
    window.setTimeout(() => citationCopied = false, 2000);
  }

  const chapters = [
    {
      name: 'Preliminaries and finite fields',
      summary: 'Introduces the algebraic language used throughout the thesis, then constructs finite fields and their algebraic closures.'
    },
    {
      name: 'Number fields and ideals',
      summary: 'Develops number rings, ideals, norms, ramification and class groups, with one example connecting the concepts.'
    },
    {
      name: 'Elliptic curves',
      summary: 'Builds the group law geometrically and algebraically before studying elliptic curves over finite fields.'
    },
    {
      name: 'Morphisms of elliptic curves',
      summary: 'Explains isogenies, isomorphism classes, Montgomery curves, endomorphism rings, ideal actions and isogeny graphs.'
    },
    {
      name: 'Encryption schemes',
      summary: 'Moves from Diffie–Hellman key exchange to CSIDH, then analyses a CSIDH variant and its security.'
    }
  ];

  const schema = {
    '@context': 'https://schema.org',
    '@type': 'ScholarlyArticle',
    headline: 'Decoding CSIDH: A Guide to Isogeny-Based Cryptography',
    name: 'Decoding CSIDH: A Guide to Isogeny-Based Cryptography',
    abstract,
    author: {
      '@type': 'Person',
      name: 'L. J. Perrenet',
      url: 'https://jorisperrenet.com/'
    },
    publisher: {
      '@type': 'CollegeOrUniversity',
      name: 'Leiden University',
      department: 'Mathematical Institute'
    },
    datePublished: '2024-07-22',
    url: canonical,
    sameAs: 'https://github.com/jorisperrenet/MasterThesis',
    inLanguage: 'en',
    associatedMedia: {
      '@type': 'MediaObject',
      contentUrl: pdfUrl,
      encodingFormat: 'application/pdf'
    },
    keywords: [
      'CSIDH',
      'isogeny-based cryptography',
      'post-quantum cryptography',
      'elliptic curves',
      'ideal class group',
      'SageMath'
    ]
  };
</script>

<SiteMeta
  {title}
  {description}
  {canonical}
  type="article"
  image="https://jorisperrenet.com/screenshots/master-thesis.png"
  imageWidth={548}
  imageHeight={543}
/>

<svelte:head>
  <meta name="author" content="L. J. Perrenet" />
  <meta name="keywords" content="CSIDH, isogeny-based cryptography, post-quantum cryptography, elliptic curves, ideal class group, supersingular curves, SageMath" />
  <meta name="citation_title" content="Decoding CSIDH: A Guide to Isogeny-Based Cryptography" />
  <meta name="citation_author" content="L. J. Perrenet" />
  <meta name="citation_publication_date" content="2024/07/22" />
  <meta name="citation_dissertation_institution" content="Leiden University" />
  <meta name="citation_pdf_url" content={pdfUrl} />
  <meta property="article:published_time" content="2024-07-22" />
  {@html `<script type="application/ld+json">${JSON.stringify(schema).replace(/</g, '\\u003c')}</script>`}
</svelte:head>

<article class="mx-auto max-w-5xl px-4 py-8 sm:px-6 sm:py-12">
  <header class="grid items-center gap-7 rounded-3xl border border-blue-100 bg-gradient-to-br from-white to-blue-50 p-5 shadow-sm dark:border-gray-700 dark:from-gray-800 dark:to-gray-900 sm:grid-cols-[220px_1fr] sm:p-8">
    <img
      src={`${base}/master-thesis.png`}
      width="548"
      height="543"
      alt="Preview of Decoding CSIDH, showing elliptic-curve constructions from the thesis"
      class="mx-auto aspect-square w-full max-w-52 rounded-xl border border-gray-200 object-cover shadow-xl dark:border-gray-600"
    />
    <div>
      <p class="m-0 text-sm font-bold uppercase tracking-[0.18em] text-blue-700 dark:text-blue-300">Master thesis · Leiden University</p>
      <h1 class="mt-3 text-3xl font-extrabold tracking-tight text-gray-950 dark:text-white sm:text-5xl">Decoding CSIDH</h1>
      <p class="mt-3 text-lg text-gray-600 dark:text-gray-300 sm:text-xl">A Guide to Isogeny-Based Cryptography</p>
      <p class="mt-4 text-sm leading-6 text-gray-600 dark:text-gray-300">L. J. Perrenet · Mathematical Institute · July 22, 2024 · Supervised by Dr. S. A. Arpin</p>
      <div class="mt-6 flex flex-wrap gap-3">
        <a class="rounded-xl bg-blue-600 px-5 py-3 font-bold text-white shadow-sm transition hover:bg-blue-700" href={`${base}/decoding-csidh.pdf`}>Read thesis PDF</a>
        <a class="rounded-xl border border-gray-300 bg-white px-5 py-3 font-bold text-gray-800 transition hover:border-blue-400 hover:text-blue-700 dark:border-gray-600 dark:bg-gray-800 dark:text-gray-100 dark:hover:border-blue-400 dark:hover:text-blue-300" href={`${base}/decoding-csidh.pdf`} download>Download PDF</a>
        <a class="rounded-xl px-3 py-3 font-semibold text-blue-700 hover:underline dark:text-blue-300" href="#citation">Cite this thesis</a>
      </div>
    </div>
  </header>

  <section aria-labelledby="abstract" class="mx-auto mt-12 max-w-3xl">
    <h2 id="abstract" class="text-2xl font-bold tracking-tight text-gray-950 dark:text-white">Abstract</h2>
    <p class="mt-4 text-base leading-8 text-gray-700 dark:text-gray-300">{abstract}</p>
  </section>

  <section aria-labelledby="guide" class="mt-12">
    <div class="mx-auto max-w-3xl">
      <p class="text-sm font-bold uppercase tracking-[0.16em] text-blue-700 dark:text-blue-300">What you will learn</p>
      <h2 id="guide" class="mt-2 text-2xl font-bold tracking-tight text-gray-950 dark:text-white">Contents and chapter guide</h2>
      <p class="mt-3 leading-7 text-gray-600 dark:text-gray-300">The thesis is organised into five main chapters, preceded by prerequisites and followed by the SageMath implementation in Appendix A.</p>
    </div>
    <ol class="mx-auto mt-6 grid max-w-3xl gap-4">
      {#each chapters as chapter, index}
        <li class="grid grid-cols-[2.5rem_1fr] gap-3 rounded-2xl border border-gray-200 bg-white p-4 dark:border-gray-700 dark:bg-gray-800">
          <span class="flex h-9 w-9 items-center justify-center rounded-full bg-blue-100 font-bold text-blue-800 dark:bg-blue-950 dark:text-blue-200">{index + 1}</span>
          <div>
            <h3 class="font-bold text-gray-950 dark:text-white">{chapter.name}</h3>
            <p class="mt-1 leading-7 text-gray-600 dark:text-gray-300">{chapter.summary}</p>
          </div>
        </li>
      {/each}
    </ol>
  </section>

  <section aria-labelledby="highlights" class="mx-auto mt-12 max-w-3xl">
    <h2 id="highlights" class="text-2xl font-bold tracking-tight text-gray-950 dark:text-white">Worked-example highlights</h2>
    <div class="mt-5 grid gap-4 sm:grid-cols-3">
      <div class="rounded-2xl bg-blue-50 p-4 dark:bg-blue-950/40">
        <h3 class="font-bold text-blue-950 dark:text-blue-100">Elliptic-curve arithmetic</h3>
        <p class="mt-2 text-sm leading-6 text-blue-900/80 dark:text-blue-200">Geometric constructions and explicit formulas make addition and doubling inspectable.</p>
      </div>
      <div class="rounded-2xl bg-violet-50 p-4 dark:bg-violet-950/40">
        <h3 class="font-bold text-violet-950 dark:text-violet-100">Ideal actions</h3>
        <p class="mt-2 text-sm leading-6 text-violet-900/80 dark:text-violet-200">Concrete examples connect ideal classes to walks through isogeny graphs.</p>
      </div>
      <div class="rounded-2xl bg-emerald-50 p-4 dark:bg-emerald-950/40">
        <h3 class="font-bold text-emerald-950 dark:text-emerald-100">CSIDH in SageMath</h3>
        <p class="mt-2 text-sm leading-6 text-emerald-900/80 dark:text-emerald-200">Appendix A translates the mathematical construction into a readable implementation.</p>
      </div>
    </div>
  </section>

  <section id="citation" aria-labelledby="citation-title" class="mx-auto mt-12 max-w-3xl scroll-mt-20 rounded-2xl border border-gray-200 bg-white p-5 dark:border-gray-700 dark:bg-gray-800 sm:p-6">
    <div class="flex items-center justify-between gap-4">
      <h2 id="citation-title" class="text-2xl font-bold tracking-tight text-gray-950 dark:text-white">BibTeX citation</h2>
      <button
        type="button"
        onclick={copyCitation}
        class="shrink-0 rounded-lg bg-gray-100 px-3 py-2 text-sm font-bold text-gray-800 transition hover:bg-gray-200 dark:bg-gray-700 dark:text-gray-100 dark:hover:bg-gray-600"
        aria-live="polite"
      >{citationCopied ? 'Copied!' : 'Copy BibTeX'}</button>
    </div>
    <pre class="mt-4 overflow-x-auto rounded-xl bg-gray-950 p-4 text-sm leading-6 text-gray-100"><code>{bibtex}</code></pre>
    <div class="mt-4 flex flex-wrap gap-x-5 gap-y-2 text-sm font-semibold">
      <a class="text-blue-700 hover:underline dark:text-blue-300" href={pdfUrl}>Open the thesis</a>
      <a class="text-blue-700 hover:underline dark:text-blue-300" href="https://github.com/jorisperrenet/MasterThesis">LaTeX and SageMath source</a>
    </div>
  </section>
</article>
